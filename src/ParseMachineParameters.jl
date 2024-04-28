using Base.Iterators
using Random
using LightXML


function dicts_to_LookUpTable(dicts, T=Float32)
    max_n_entries = maximum([length(d["values"]) for d in dicts])
    energies = Vector{T}(undef, length(dicts))
    x0 = Vector{T}(undef, length(dicts))
    dx = Vector{T}(undef, length(dicts))
    table = Matrix{T}(undef, length(dicts), max_n_entries)
    fill!(table, NaN)

    for (i, d) in enumerate(dicts)
        energies[i] = parse(T, d["en"])
        x0[i] = parse(T, d["x0"])
        dx[i] = parse(T, d["dx"])
        table[i, 1:length(d["values"])] .= d["values"]
    end

    permutation = sortperm(energies)
    energies = energies[permutation]
    x0 = x0[permutation]
    dx = dx[permutation]
    table = table[permutation, :]

    @assert issorted(energies)

    # Precalculate where the lookup table for each energy ends.
    n = length(energies)
    indices_end = Vector{Int32}(undef, n)
    for i in 1:n
        ind = findfirst(isnan, table[i, :])
        if isnothing(ind)
            indices_end[i] = size(table, 2)
        else
            indices_end[i] = ind
        end
    end

    # Convert to keV.
    return Juliana.LookUpTable(energies * 1000, x0, dx, table, indices_end)
end


function parse_beamdata(path, T=Float32)
    all_dose_dicts = []
    all_mcs_dicts = []
    all_let_dicts = []

    doc = LightXML.parse_file(path);
    dat = LightXML.root(doc)
    tables = child_elements(dat);

    for table in tables
        dict = attributes_dict(table)
        dict = convert(Dict{String, Any}, dict)
        dict["values"] = Vector{T}(undef, 0)
        lut = Iterators.only([c for c in child_elements(table)])
        for entry in child_elements(lut)
            @assert LightXML.name(entry) == "entry"
            push!(dict["values"], parse(T, content(entry)))
        end

        if LightXML.name(table) == "dose"
            push!(all_dose_dicts, dict)
        elseif LightXML.name(table) == "mcs"
            push!(all_mcs_dicts, dict)
        elseif LightXML.name(table) == "let"
            push!(all_let_dicts, dict)
        else
            @error "Unknown element name: $(LightXML.name(table))"
        end
    end
    
    depth_dose_curves = dicts_to_LookUpTable(all_dose_dicts);
    sigma_mcs_curves = dicts_to_LookUpTable(all_mcs_dicts);
    let_curves = dicts_to_LookUpTable(all_let_dicts);
    
    return depth_dose_curves, sigma_mcs_curves, let_curves
end


function parse_symmetric_phasespace(path, T=Float32)
    doc = LightXML.parse_file(path);
    dat = LightXML.root(doc)
    all_params = collect(child_elements(dat));

    energies = Vector{T}(undef, length(all_params))
    a0 = Vector{T}(undef, length(all_params))
    a1 = Vector{T}(undef, length(all_params))
    a2 = Vector{T}(undef, length(all_params))

    for (i, params) in enumerate(all_params)
        @assert LightXML.name(params) == "params"
        dict = attributes_dict(params)
        energies[i] = parse(T, dict["en"]) * 1000 # convert to keV
        a0[i] = parse(T, dict["a0"])
        a1[i] = parse(T, dict["a1"])
        a2[i] = parse(T, dict["a2"])
    end
    
    permutation = sortperm(energies)
    energies = energies[permutation]
    a0 = a0[permutation]
    a1 = a1[permutation]
    a2 = a2[permutation]
    
    @assert issorted(energies)

    return Juliana.PhaseSpace(energies, a0, a1, a2)
end


function parse_asymmetric_phasespace(path, T=Float32)
    doc = LightXML.parse_file(path);
    dat = LightXML.root(doc)
    all_params = collect(child_elements(dat));

    energies = Vector{T}(undef, length(all_params))
    a0t = Vector{T}(undef, length(all_params))
    a1t = Vector{T}(undef, length(all_params))
    a2t = Vector{T}(undef, length(all_params))
    a0u = Vector{T}(undef, length(all_params))
    a1u = Vector{T}(undef, length(all_params))
    a2u = Vector{T}(undef, length(all_params))

    for (i, params) in enumerate(all_params)
        @assert LightXML.name(params) == "params"
        dict = attributes_dict(params)
        energies[i] = parse(T, dict["en"]) * 1000 # convert to keV
        a0t[i] = parse(T, dict["a0t"])
        a1t[i] = parse(T, dict["a1t"])
        a2t[i] = parse(T, dict["a2t"])
        a0u[i] = parse(T, dict["a0u"])
        a1u[i] = parse(T, dict["a1u"])
        a2u[i] = parse(T, dict["a2u"])
    end
    
    permutation = sortperm(energies)
    energies = energies[permutation]
    a0t = a0t[permutation]
    a1t = a1t[permutation]
    a2t = a2t[permutation]
    a0u = a0u[permutation]
    a1u = a1u[permutation]
    a2u = a2u[permutation]
    
    @assert issorted(energies)

    return Juliana.PhaseSpace(energies, a0t, a1t, a2t, a0u, a1u, a2u, false)
end


function build_phase_space(depth_dose_lut,
                           beamline_phasespace,
                           absorber_phasespace,
                           nozzle_phasespace,
                           nozzle_extraction,
                           n_preabsorbers,
                           T=Float32)
    @assert depth_dose_lut.energies == beamline_phasespace.energies
    @assert depth_dose_lut.energies == absorber_phasespace.energies
    @assert issorted(depth_dose_lut.energies)
    @assert issorted(nozzle_phasespace.energies)
    @assert nozzle_phasespace.is_symmetric

    maxDistance = 46.0f0
    nozzleIsoDistance = maxDistance - nozzle_extraction;

    energies = depth_dose_lut.energies
    n_energies = length(energies)

    phase_space_a0t = Vector{T}(undef, n_energies)
    phase_space_a1t = Vector{T}(undef, n_energies)
    phase_space_a2t = Vector{T}(undef, n_energies)
    phase_space_a0u = Vector{T}(undef, n_energies)
    phase_space_a1u = Vector{T}(undef, n_energies)
    phase_space_a2u = Vector{T}(undef, n_energies)

    for (i, E) in enumerate(energies)
        nozzle_index = argmin(abs.(E .- nozzle_phasespace.energies))

        # We could as well take the values in u direction because the phase space is symmetric.
        a0 = nozzle_phasespace.a0t[nozzle_index]
        a1 = nozzle_phasespace.a0t[nozzle_index] * nozzleIsoDistance   +     nozzle_phasespace.a1t[nozzle_index]
        a2 = nozzle_phasespace.a0t[nozzle_index] * nozzleIsoDistance^2 + 2 * nozzle_phasespace.a1t[nozzle_index] * nozzleIsoDistance + nozzle_phasespace.a2t[nozzle_index]

        sigma_index = argmin(abs.(E .- beamline_phasespace.energies))
        a0t = beamline_phasespace.a0t[sigma_index] + a0;
        a1t = beamline_phasespace.a1t[sigma_index] + a1;
        a2t = beamline_phasespace.a2t[sigma_index] + a2;
        a0u = beamline_phasespace.a0u[sigma_index] + a0;
        a1u = beamline_phasespace.a1u[sigma_index] + a1;
        a2u = beamline_phasespace.a2u[sigma_index] + a2;

        if n_preabsorbers > 0
            @assert absorber_phasespace.is_symmetric
            absorber_index = argmin(abs.(E .- absorber_phasespace.energies))
            a0_preabsorber = absorber_phasespace.a0t[absorber_index]
            a1_preabsorber = absorber_phasespace.a0t[absorber_index] * nozzleIsoDistance   +     absorber_phasespace.a1t[absorber_index]
            a2_preabsorber = absorber_phasespace.a0t[absorber_index] * nozzleIsoDistance^2 + 2 * absorber_phasespace.a1t[absorber_index] * nozzleIsoDistance + absorber_phasespace.a2t[absorber_index]

            a0t += a0_preabsorber
            a1t += a1_preabsorber
            a2t += a2_preabsorber
            a0u += a0_preabsorber
            a1u += a1_preabsorber
            a2u += a2_preabsorber
        end

        phase_space_a0t[i] = a0t
        phase_space_a1t[i] = a1t
        phase_space_a2t[i] = a2t
        phase_space_a0u[i] = a0u
        phase_space_a1u[i] = a1u
        phase_space_a2u[i] = a2u
    end
    
    return Juliana.PhaseSpace(
        convert.(T, collect(energies)),
        phase_space_a0t,
        phase_space_a1t,
        phase_space_a2t,
        phase_space_a0u,
        phase_space_a1u,
        phase_space_a2u,
        false,
    )
end


function load_machine_parameters(directory, nozzle_extraction)
    path = "$(directory)/beamdata.xml"
    depth_dose_curves, sigma_mcs_curves, let_curves = parse_beamdata(path)
    
    path = "$(directory)/nozzle.xml"
    nozzle_phasespace = parse_symmetric_phasespace(path)

    path = "$(directory)/absorber.xml"
    absorber_phasespace = parse_symmetric_phasespace(path)

    path = "$(directory)/beamline.xml"
    beamline_phasespace = parse_asymmetric_phasespace(path)
    
    phase_space_no_preabsorber = build_phase_space(
        depth_dose_curves,
        beamline_phasespace,
        absorber_phasespace,
        nozzle_phasespace,
        nozzle_extraction,
        0,
    )
    phase_space_with_preabsorber = build_phase_space(
        depth_dose_curves,
        beamline_phasespace,
        absorber_phasespace,
        nozzle_phasespace,
        nozzle_extraction,
        1,
    )
    
    return depth_dose_curves, sigma_mcs_curves, phase_space_no_preabsorber, phase_space_with_preabsorber
end
