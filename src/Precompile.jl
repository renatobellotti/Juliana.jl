# This file aims to precompile some often used function to speed up loading
# the execution. Details can be found in the documentation of the
# SnoopCompile.jl package:
# https://timholy.github.io/SnoopCompile.jl/dev/snoopi_deep_parcel/#precompilation
using CUDA


# if we're precompiling the package
if ccall(:jl_generating_output, Cint, ()) == 1
    # let
    #     T = Float32

    #     vector_types = (
    #         Vector{T},
    #         # CuVector{T},
    #     )

    #     matrix_types = (
    #         Matrix{T},
    #         CuMatrix{T},
    #     )

    #     for (VT, MT) in zip(vector_types, matrix_types)
    #         v = VT(undef, 1)
    #         Dij = MT(undef, 1, 1)
    #         fill!(v, 1)
    #         fill!(Dij, 1)

    #         presc = Juliana.Prescriptions{T}(
    #             Vector{Tuple{String, T}}(undef, 0),
    #             Vector{Juliana.Constraint{T}}(undef, 0),
    #         )
    #         structures = Dict{String, VT}()

    #         config = Juliana.OptimisationConfiguration(
    #             1f0, # normalisationDose
    #             v,   # normalisationStructureMask
    #             v,   # ct
    #             v,   # idealDose
    #             v,   # importance
    #             v,   # minimumDose
    #             Dij, # Dij,
    #             structures,
    #             presc,
    #         );

    #         dose = VT(undef, 1)
    #         safety_margin = 0f0

    #         subloss_weights = Dict{String, T}()
            
    #         Juliana.build_loss_parts(dose, config, safety_margin)
    #         Juliana.dose_loss(dose, config, subloss_weights)
    #     end
    # end
end