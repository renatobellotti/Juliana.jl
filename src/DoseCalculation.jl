# The source code in this file is a port from the Fiona treatment planning
# system used at the Center for Proton Therapy at the Paul Scherrer Institut.
# It was ported by Bastien Golomer. The resulting code was adjusted by
# Renato Bellotti to better match the package.

using Atomix
using CUDA
using CUDAKernels
using KernelAbstractions


# A constant in Fiona, so we just keep it  as it is.
const PREABSORBER_WED =  4.173f0
const MIN_DIJ_CONTRIBUTION = convert(Float32, 1.e-8)
const HU_AIR = -900


"""
    computeDose(stuFieldCenter, wed, energy, tSpot, uSpot, n_preabsorbers, ddc, sigma_mcs_lut, phase_space_no_preabsorber, phase_space_with_preabsorber, preabsorberWed) 

This function computes the dose at a selected optimisation point from a
selected spot. 
"""
function computeDose(stuFieldCenter,
                     # Characterisation of the optimisation point.
                     wed,
                     # Spot variables.
                     energy,       # energy to be used for the depth-dose lookup
                     energy_sigma, # energy used for the beam size lookup
                     tSpot,
                     uSpot,
                     n_preabsorbers,
                     # Characteristics of the beam physics.
                     ddc,
                     sigma_mcs_lut,
                     phase_space_no_preabsorber,
                     phase_space_with_preabsorber,
                     preabsorberWed) 
    # Constants from Fiona.
    CUTOFF_FACTOR_SQUARED = 25.f0
    # End of constants.

    s, t, u = stuFieldCenter

    if n_preabsorbers == 0
        a0t, a1t, a2t, a0u, a1u, a2u = getPhaseSpaceParameters(phase_space_no_preabsorber, energy)
    elseif n_preabsorbers == 1
        a0t, a1t, a2t, a0u, a1u, a2u = getPhaseSpaceParameters(phase_space_with_preabsorber, energy)
    end

    aIPSt = a2t + 2 * a1t * s + a0t * s^2;
    aIPSu = a2u + 2 * a1u * s + a0u * s^2;

    tDistance = t - tSpot;
    uDistance = u - uSpot;

    result = 0.f0
    if ((tDistance^2 + uDistance^2) < (CUTOFF_FACTOR_SQUARED * aIPSt / 2.f0)) && (1.e-10 < wed)
        longitudinal_dose = interpolate_lut(
            wed + n_preabsorbers * preabsorberWed,
            energy,
            ddc,
        )

        lateral_dose = lateral_dose_curve(
            sigma_mcs_lut,
            energy_sigma,
            wed, # Not a typo: the MCS from the preabsorber is already in the energy
            tDistance,
            uDistance,
            aIPSt,
            aIPSu,
        )

        result = longitudinal_dose * lateral_dose
    end

    return result
end


@kernel function DijInterpolationKernel(Dij,
                                        @Const(wed),
                                        # Variables characterising the spots.
                                        @Const(field_center),
                                        @Const(gantry_angle),
                                        @Const(couch_angle),
                                        @Const(spots_t),
                                        @Const(spots_u),
                                        @Const(spots_energies),
                                        @Const(spots_absorbers),
                                        # Technical parameters.
                                        @Const(optimisation_points),
                                        # Machine characteristics.
                                        @Const(ddc),
                                        @Const(sigma_mcs_lut),
                                        @Const(phase_space_no_preabsorber),
                                        @Const(phase_space_with_preabsorber),
                                        @Const(preabsorberWed))
    optIndex = @index(Global, Cartesian)[1]
    spotIndex = @index(Global, Cartesian)[2]

    n_preabsorbers = spots_absorbers[spotIndex]
    t = spots_t[spotIndex]
    u = spots_u[spotIndex]
    energy = spots_energies[spotIndex]
    energy_sigma = energy_after_preabsorber(energy, n_preabsorbers * preabsorberWed)

    gantry_angle = gantry_angle
    couch_angle = couch_angle

    optimPointX = optimisation_points[optIndex, 1]
    optimPointY = optimisation_points[optIndex, 2]
    optimPointZ = optimisation_points[optIndex, 3]

    stuFieldCenter = XyzToStu(
        optimPointX,
        optimPointY,
        optimPointZ,
        field_center, 
        gantry_angle,
        couch_angle,
    )

    d = computeDose(
        stuFieldCenter,
        wed[optIndex],
        energy,
        energy_sigma,
        t,
        u,
        n_preabsorbers,
        ddc,
        sigma_mcs_lut,
        phase_space_no_preabsorber,
        phase_space_with_preabsorber,
        preabsorberWed,
    )
    if d < MIN_DIJ_CONTRIBUTION
        d = 0
    end
    Dij[optIndex, spotIndex] = d
end


dij_kernel = DijInterpolationKernel(CUDADevice(), 32)


function runDijKernel(Dij, 
                      wed,
                      field_center,
                      gantry_angle,
                      couch_angle,
                      spots_t,
                      spots_u,
                      spots_energies,
                      spots_absorbers,
                      optimisation_points,
                      ddc,
                      sigma_mcs_lut,
                      phase_space_no_preabsorber,
                      phase_space_with_preabsorber;
                      preabsorberWed=PREABSORBER_WED)
    event = dij_kernel(
        Dij, 
        wed,
        field_center,
        gantry_angle,
        couch_angle,
        spots_t,
        spots_u,
        spots_energies,
        spots_absorbers,
        optimisation_points,
        ddc,
        sigma_mcs_lut,
        phase_space_no_preabsorber,
        phase_space_with_preabsorber,
        preabsorberWed;
        ndrange=size(Dij),
    )
    wait(event);
end


@kernel function DoseCalculationKernel(dose,
                                       @Const(wed),
                                       # Variables characterising the spots.
                                       @Const(field_center),
                                       @Const(gantry_angle),
                                       @Const(couch_angle),
                                       @Const(spots_t),
                                       @Const(spots_u),
                                       @Const(spots_energies),
                                       @Const(spots_absorbers),
                                       @Const(spots_weights),
                                       # Technical parameters.
                                       @Const(optimisation_points),
                                       # Machine characteristics.
                                       @Const(ddc),
                                       @Const(sigma_mcs_lut),
                                       @Const(phase_space_no_preabsorber),
                                       @Const(phase_space_with_preabsorber),
                                       @Const(preabsorberWed))
    optIndex = @index(Global, Cartesian)[1]

    gantry_angle = gantry_angle
    couch_angle = couch_angle

    optimPointX = optimisation_points[optIndex, 1]
    optimPointY = optimisation_points[optIndex, 2]
    optimPointZ = optimisation_points[optIndex, 3]

    stuFieldCenter = XyzToStu(
        optimPointX,
        optimPointY,
        optimPointZ,
        field_center, 
        gantry_angle,
        couch_angle,
    )

    for spotIndex in 1:length(spots_weights)
        n_preabsorbers = spots_absorbers[spotIndex]
        t = spots_t[spotIndex]
        u = spots_u[spotIndex]
        energy = spots_energies[spotIndex]
        w = spots_weights[spotIndex]
        energy_sigma = energy_after_preabsorber(energy, n_preabsorbers * preabsorberWed)

        d = computeDose(
            stuFieldCenter,
            wed[optIndex],
            energy,
            energy_sigma,
            t,
            u,
            n_preabsorbers,
            ddc,
            sigma_mcs_lut,
            phase_space_no_preabsorber,
            phase_space_with_preabsorber,
            preabsorberWed,
        )
        dose[optIndex] += w * d
    end
end


dose_kernel = DoseCalculationKernel(CUDADevice(), 32)


function getPhaseSpaceParameters(phase_space, energy)
    energyIndex = getEnergyIndex(energy, phase_space.energies)

    return phase_space.a0t[energyIndex], phase_space.a1t[energyIndex], phase_space.a2t[energyIndex],
           phase_space.a0u[energyIndex], phase_space.a1u[energyIndex], phase_space.a2u[energyIndex]
end


"""
    getNonNaNLength(vector::CuArray{Float32, 2},index::Int32)

Computes the length of the input array up until the first NaN value.

# Arguments
- `vector::CuArray{Float32, 2}`: datastructure containing all the arrays of interest. We select one array using next argument
- `index::Int32`: index for the desired array in the previously  definied datastructure

# Returns
- `nonNaNLength::Int32`: length until a NaN value

"""
function getNonNaNLength(vector, index::Integer)
    nonNaNLength = 0
    while(!isnan(vector[index,nonNaNLength+1]))
        nonNaNLength += 1
    end
    return nonNaNLength
end


"""
    function energy_after_preabsorber(energykeV::Float32, absorberThickness::Float32)

Calculate the energy of the beam after the preabsorber (if any).

The code relies on Bortfeld 1997, https://doi.org/10.1118/1.598116.
"""
function energy_after_preabsorber(energykeV::Float32, absorberThickness::Float32)
    ALPHA = 0.00220f0
    P = 1.77f0

    if (absorberThickness < 1e-6)
        return energykeV;
    else
        range = ALPHA * (energykeV / 1000)^P
        range = clamp(range, absorberThickness, Inf32)
        energy_keV = 1000 * ALPHA^(-1 / P) * (range - absorberThickness)^(1 / P)
        return energy_keV
    end
end;


"""
    getEnergyIndex(energykeV::Float32, energies::CuArray{Float32, 1})

matches an input energy with a reference table and gets the index of the matching row in the reference table.

# Arguments
- `energykeV::Float32`: enregy value to match
- `energies::CuArray{Float32, 1}`: table of energies 

# Returns
- `index::Int32`: index of the matching row in the reference table

"""
function getEnergyIndex(energykeV::Float32, energies)::Int32
    # OK, careful units of energies
    THRESHOLD = 0.1f0
    index = 0
    minDifference = Inf32
    
    # here we select the right index, and it is consistent with 1 based indexing
    for i in 1:length(energies)
        difference = abs(energies[i] - energykeV)
        if (difference < THRESHOLD) 
            return Int32(i)
        elseif (difference < minDifference) 
            minDifference = difference;
            index = i
        end
    end
    return index
end


function interpolate_lut(depth::Real,
                         energy::Real,
                         ddc::LookUpTable)
    energyIndex = getEnergyIndex(energy, ddc.energies)
    
    depth0 = ddc.x0[energyIndex]
    depthSpacing = ddc.dx[energyIndex]
    
    # +1 since Julia is 1 based indexed, has to be floor to replicate fiona
    index0 = Int(floor((depth - depth0) / depthSpacing)) + 1

    nonNaNLength = getNonNaNLength(ddc.table, energyIndex)

    if (index0 >= nonNaNLength)
        # Extrapolate using a constant if the depth larger than
        # the largest value in the lookup table (LUT).
        return ddc.table[energyIndex, nonNaNLength]
    elseif index0 <= 0
        # Extrapolate using a constant if the depth smaller than
        # the smallest value in the lookup table (LUT).
        return ddc.table[energyIndex, 1];  
        
    else
        # Interpolate dose as a function of depth linearly.
        # -1: Julia indices start at 1.
        depthFactor = (depth - (depth0 + (index0-1) * depthSpacing)) / depthSpacing;
        dose0 = ddc.table[energyIndex, index0]
        dose1 = ddc.table[energyIndex, index0+1];
        return dose0 + (dose1 - dose0) * depthFactor;
    end
    return 0f0
end


function lateral_dose_curve(sigma_mcs_lut,
                            energy,
                            wed,
                            tDistance,
                            uDistance,
                            aIPSt,
                            aIPSu)
    sigmaMcs = interpolate_lut(wed, energy, sigma_mcs_lut)
    aMcs = 2 * sigmaMcs^2;

    at =  aMcs + aIPSt;
    au =  aMcs + aIPSu;

    return exp(-tDistance^2 / at -  uDistance^2 / au) / max((π * sqrt(at * au)), 0.000001f0)
end


################################
# Dose calculation on a grid.
################################
@kernel function dose_calculation_kernel(d_points_stu,
                                         d_spots,
                                         d_wed_linear,
                                         d_dose_linear,
                                         d_depth_dose_curves,
                                         d_sigma_mcs_curves,
                                         d_phase_space_no_preabsorber,
                                         d_phase_space_with_preabsorber)
    index = @index(Global, Cartesian)
    i = index[1]
    j = index[2]

    point_stu = (
        d_points_stu[1, i],
        d_points_stu[2, i],
        d_points_stu[3, i],
    )
    wed = d_wed_linear[i]
    spot = d_spots[j]

    n_preabsorbers = spot.numberOfAbsorbers
    E = spot.energykeV
    E_sigma = Juliana.energy_after_preabsorber(E, n_preabsorbers * Juliana.PREABSORBER_WED)

    spot_dose = Juliana.computeDose(
        point_stu,
        wed,
        E,
        E_sigma,
        spot.t,
        spot.u,
        n_preabsorbers,
        d_depth_dose_curves,
        d_sigma_mcs_curves,
        d_phase_space_no_preabsorber,
        d_phase_space_with_preabsorber,
        Juliana.PREABSORBER_WED,
    )

    Atomix.@atomic d_dose_linear[i] += spot.weight * spot_dose
end


dose_calculation_gpu = dose_calculation_kernel(CUDADevice(), 32);


"""
    calculate_dose(field, points, wed_linear, tps::Juliana.JulianaTps)

Calcuate the dose for the given field at the given points with precalculated WED.
"""
function calculate_dose(field, points, wed_linear, tps::Juliana.JulianaTps)
    n_points = size(points, 2)
    n_spots = length(field.spots)

    # Calculate points of interest in STU coordinates.
    field_center = (
        field.fieldCenter["x"],
        field.fieldCenter["y"],
        field.fieldCenter["z"],
    )

    points_stu = Juliana.XyzToStu(
        points,
        field_center,
        field.gantryAngle,
        field.couchAngle,
    )

    # Allocate memory for the output.
    dose_linear = Vector{Float32}(undef, n_points)
    fill!(dose_linear, 0)

    # Move data to the GPU.
    d_points_stu = cu(points_stu)
    d_spots = cu(field.spots)
    d_wed_linear = cu(wed_linear)
    d_dose_linear = cu(dose_linear)
    d_tps = Juliana.to_gpu(tps)

    # Call the dose calculation on the GPU.
    event = dose_calculation_gpu(
        d_points_stu,
        d_spots,
        d_wed_linear,
        d_dose_linear,
        d_tps.d_depth_dose_curves,
        d_tps.d_sigma_mcs_curves,
        d_tps.d_phase_space_no_preabsorber,
        d_tps.d_phase_space_with_preabsorber,
        ndrange=(n_points, n_spots),
    )
    wait(event)

    # Move back to the CPU.
    dose_linear = collect(d_dose_linear)

    return dose_linear
end


function calculate_dose(tps::Juliana.JulianaTps,
                        resolution::Real,
                        ct::Juliana.ScalarGrid,
                        plan::Juliana.TreatmentPlan)
    dose_calc_grid = Juliana.get_fiona_dose_calc_grid(
        ct.grid,
        resolution,
    )
    points, indices = Juliana.grid_to_points_and_indices(dose_calc_grid)

    ct_dose_calc_grid = Juliana.interpolate_linearly(ct, dose_calc_grid)
    point_not_air = Juliana.volume_at_indices(ct_dose_calc_grid.data, indices) .> HU_AIR

    points = points[:, point_not_air]
    indices = indices[:, point_not_air]

    densities = Juliana.ScalarGrid(
        tps.convert_to_sp.(ct.data),
        ct.grid,
    )

    field_doses = []
    for field in plan.fields
        push!(field_doses, calculate_dose(
            tps,
            resolution,
            densities,
            points,
            field,
        ))
    end
    dose = .+(field_doses...)

    dose_cube = Juliana.ScalarGrid(
        convert.(Float32, Juliana.flat_vector_to_cube(dose_calc_grid, indices, dose)),
        dose_calc_grid,
    )
    not_air_cube = Juliana.ScalarGrid(
        # Select only the not-air voxels to match with the indices array.
        convert.(Float32, Juliana.flat_vector_to_cube(dose_calc_grid, indices, point_not_air[point_not_air])),
        dose_calc_grid,
    )
    dose_cube = interpolate_linearly(dose_cube, not_air_cube, ct.grid)

    return dose_cube
end


function calculate_dose(tps::Juliana.JulianaTps,
                        resolution::Real,
                        densities,
                        points,
                        field::Juliana.FieldDefinition)
    wed = vec(Juliana.calculate_wed(
        densities,
        [field.gantryAngle],
        [field.couchAngle],
        points,
    ))
    field_dose = calculate_dose(
        field,
        points,
        wed,
        tps,
    )
    return field_dose
end
