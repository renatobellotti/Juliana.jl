# The NTCP models in this file are a direct port of the Excel sheet that is
# used at the center for proton therapy at PSI (version 3.0).
# They are based on the following paper:
# https://doi.org/10.1016/j.radonc.2011.05.010

function ntcp(S)
    return 1 / (1 + exp(-S))
end


function xerostomia_parameters(grade, baseline_xerostomia_text, setting)
    @assert baseline_xerostomia_text in ["none", "mild", "moderate", "severe"]
    @assert grade in [2, 3]
    @assert setting in ["primary", "post-operative"]

    baseline = 0
    offset = 0
    ipsilateral_parotids = 0
    submandibulars = 0
    baseline_xerostomia = 0

    if setting == "primary"
        if grade == 2 # grade ≥ 2
            if baseline_xerostomia_text == "none"
                baseline = 0
            elseif baseline_xerostomia_text == "mild"
                baseline = 0.495
            else
                baseline = 1.207
            end
            offset = -2.2951
            ipsilateral_parotids = 0.0996
            submandibulars = 0.0182
            baseline_xerostomia = baseline
        else # grade ≥ 3
            if baseline_xerostomia_text == "none"
                baseline = 0
            elseif baseline_xerostomia_text == "mild"
                baseline = 0.4249
            else
                baseline = 1.0361
            end
            offset = -3.7286
            ipsilateral_parotids = 0.0855
            submandibulars = 0.0156
            baseline_xerostomia = baseline
        end
    else
        if grade == 2 # grade ≥ 2
            if baseline_xerostomia_text == "none"
                baseline = 0
            elseif baseline_xerostomia_text == "mild"
                baseline = 0.1925
            else
                baseline = 0.4695
            end
            offset = -1.6824
            ipsilateral_parotids = 0.0388
            submandibulars = 0.0071
            baseline_xerostomia = baseline
        else # grade ≥ 3
            if baseline_xerostomia_text == "none"
                baseline = 0
            elseif baseline_xerostomia_text == "mild"
                baseline = 0.5234
            else
                baseline = 1.2763
            end
            offset = -4.3613
            ipsilateral_parotids = 0.1054
            submandibulars = 0.0193
            baseline_xerostomia = baseline
        end
    end
    
    return (
        offset=offset,
        ipsilateral_parotids=ipsilateral_parotids,
        submandibulars=submandibulars,
        baseline_xerostomia=baseline_xerostomia,
    )
end


function dysphagia_parameters(grade, baseline_dysphagia_text, tumor_location_text, setting)
    @assert baseline_dysphagia_text in [0, 1, 2, 3, 4, 5]
    @assert grade in [2, 3]
    @assert tumor_location_text in ["oral_cavity", "pharynx", "larynx"]
    @assert setting in ["primary", "post-operative"]

    baseline = NaN
    offset = NaN
    oral_cavity = NaN
    pcm_superior = NaN
    pcm_medius = NaN
    pcm_inferior = NaN
    baseline_dysphagia = NaN
    tumor_location = NaN

    if setting == "primary"
        if grade == 2 # grade ≥ 2
            if baseline_dysphagia_text in [0, 1]
                baseline = 0
            elseif baseline_dysphagia_text == 2
                baseline = 0.9382
            else
                baseline = 1.29
            end
            
            if tumor_location_text == "oral_cavity"
                tumor_location = 0
            elseif tumor_location_text == "pharynx"
                tumor_location = -0.6281
            elseif tumor_location_text == "larynx"
                tumor_location = -0.7711
            else
                @error "Unknown tumor location <$tumor_location_text>. Allowed values: oral_cavity, pharynx, larynx"
            end
            
            offset = -4.0536
            oral_cavity = 0.03
            pcm_superior = 0.0236
            pcm_medius = 0.0095
            pcm_inferior = 0.0133
        else # grade ≥ 3
            if baseline_dysphagia_text in [0, 1]
                baseline = 0
            elseif baseline_dysphagia_text == 2
                baseline = 0.5738
            else
                baseline = 1.4718
            end
            
            if tumor_location_text == "oral_cavity"
                tumor_location = 0
            elseif tumor_location_text == "pharynx"
                tumor_location = 0.0387
            elseif tumor_location_text == "larynx"
                tumor_location = -0.5303
            else
                @error "Unknown tumor location <$tumor_location_text>. Allowed values: oral_cavity, pharynx, larynx"
            end
            
            offset = -7.6174
            oral_cavity = 0.0259
            pcm_superior = 0.0203
            pcm_medius = 0.0303
            pcm_inferior = 0.0341
        end
    else
        if grade == 2 # grade ≥ 2
            if baseline_dysphagia_text in [0, 1]
                baseline = 0
            elseif baseline_dysphagia_text == 2
                baseline = 0.5985
            else
                baseline = 0.8227
            end
            
            if tumor_location_text == "oral_cavity"
                tumor_location = 0
            elseif tumor_location_text == "pharynx"
                tumor_location = -0.4007
            elseif tumor_location_text == "larynx"
                tumor_location = -0.4918
            else
                @error "Unknown tumor location <$tumor_location_text>. Allowed values: oral_cavity, pharynx, larynx"
            end
            
            offset = -2.4138
            oral_cavity = 0.0192
            pcm_superior = 0.0151
            pcm_medius = 0.006
            pcm_inferior = 0.0085
        else # grade ≥ 3
            if baseline_dysphagia_text in [0, 1]
                baseline = 0
            elseif baseline_dysphagia_text == 2
                baseline = 0.1404
            else
                baseline = 0.3603
            end
            
            if tumor_location_text == "oral_cavity"
                tumor_location = 0
            elseif tumor_location_text == "pharynx"
                tumor_location = 0.0095
            elseif tumor_location_text == "larynx"
                tumor_location = -0.1298
            else
                @error "Unknown tumor location <$tumor_location_text>. Allowed values: oral_cavity, pharynx, larynx"
            end
            
            offset = -3.2594
            oral_cavity = 0.0063
            pcm_superior = 0.005
            pcm_medius = 0.0074
            pcm_inferior = 0.0084
        end
    end
    
    return (
        offset=offset,
        oral_cavity=oral_cavity,
        pcm_superior=pcm_superior,
        pcm_medius=pcm_medius,
        pcm_inferior=pcm_inferior,
        baseline_dysphagia=baseline,
        tumor_location=tumor_location,
    )
end


function xerostomia_S(mean_dose_ipsilateral_parotid,
                      mean_dose_contralateral_parotid,
                      mean_dose_submandibular_glands,
                      grade,
                      baseline_xerostomia,
                      setting)
    offset, ipsilateral_parotids, submandibulars, baseline_xerostomia = xerostomia_parameters(
        grade,
        baseline_xerostomia,
        setting,
    )

    S = offset + ipsilateral_parotids * (√(mean_dose_ipsilateral_parotid) + √(mean_dose_contralateral_parotid)) + submandibulars * mean_dose_submandibular_glands + baseline_xerostomia
    return S
end


function dysphagia_S(mean_dose_oral_cavity,
                     mean_dose_pcm_superior,
                     mean_dose_pcm_medius,
                     mean_dose_pcm_inferior,
                     grade,
                     baseline_dysphagia_text,
                     tumor_location_text,
                     setting)
    offset, oral_cavity, pcm_superior, pcm_medius, pcm_inferior, baseline_dysphagia, tumor_location = dysphagia_parameters(
        grade,
        baseline_dysphagia_text,
        tumor_location_text,
        setting,
    )
    S = offset + oral_cavity * mean_dose_oral_cavity + pcm_superior * mean_dose_pcm_superior + pcm_medius * mean_dose_pcm_medius + pcm_inferior * mean_dose_pcm_inferior + baseline_dysphagia + tumor_location
    return S
end


function xerostomia_ntcp_factory(grade, baseline_xerostomia, setting)
    @assert baseline_xerostomia in ["none", "mild", "moderate", "severe"]
    @assert grade in [2, 3]
    @assert setting in ["primary", "post-operative"]

    function xerostomia_ntcp(mean_dose_ipsilateral_parotid,
                             mean_dose_contralateral_parotid,
                             mean_dose_submandibular_glands)
        S = xerostomia_S(
            mean_dose_ipsilateral_parotid,
            mean_dose_contralateral_parotid,
            mean_dose_submandibular_glands,
            grade,
            baseline_xerostomia,
            setting,
        )
        return ntcp(S)
    end

    return xerostomia_ntcp
end


function dysphagia_ntcp_factory(grade, baseline_dysphagia, tumor_location, setting)
    @assert baseline_dysphagia in [0, 1, 2, 3, 4, 5]
    @assert grade in [2, 3]
    @assert tumor_location in ["oral_cavity", "pharynx", "larynx"]
    @assert setting in ["primary", "post-operative"]

    function dysphagia_ntcp(mean_dose_oral_cavity,
                            mean_dose_pcm_superior,
                            mean_dose_pcm_medius,
                            mean_dose_pcm_inferior)
        S = dysphagia_S(
            mean_dose_oral_cavity,
            mean_dose_pcm_superior,
            mean_dose_pcm_medius,
            mean_dose_pcm_inferior,
            grade,
            baseline_dysphagia,
            tumor_location,
            setting,
        )
        return ntcp(S)
    end

    return dysphagia_ntcp
end
