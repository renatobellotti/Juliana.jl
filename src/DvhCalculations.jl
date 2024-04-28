using Interpolations


function my_percentile(x::AbstractArray{T, 1}, percentage::T) where {T}
    x = sort(x)
    i = max(round(percentage / convert(T, 100) * length(x)), one(T))
    i = convert(Integer, i)
    return x[i]
end


function my_percentile(x::AbstractArray{T, 1}, percentage::AbstractArray{T, 1}) where {T}
    return [my_percentile(x, p) for p in percentage]
end


function dvh_d(dose, mask, v::T) where {T}
    """
    Calculates D_<v>%, i. e. the dose to the hottest v% of the voxels in the
    given structure.
    """
    # We copy to the host because quantile needs access by element,
    # which is very slow on the GPU.
    mask_array = Array(mask)
    structure_dose = reduce(vcat, Array(dose)[mask_array .== 1])
    p = my_percentile(structure_dose, v)
    return p
end


function dvh_v(dose, mask, d)
    """
    Calculates an approximation to V_<d>Gy, i. e. the volume that receives
    at least <d>Gy dose.
    """
    volume_fractions = collect(LinRange(0.f0, 100.f0, 401))
    dose_values = dvh_d(dose, mask, volume_fractions)

    Interpolations.deduplicate_knots!(dose_values)
    f = Interpolations.linear_interpolation(dose_values, volume_fractions[end:-1:1], extrapolation_bc=Interpolations.Flat())
    return f(d)
end
