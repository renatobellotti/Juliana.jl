# Important features

The code in this documentation illustrates the most interesting functions in `JulianA.jl`.

JulianA comes with an integration of the publicly available TROTS dataset. At the moment, only a single patient
is supported (with ID `"Protons1"`).

The following command will download the data from the internet and load the patient's data structure:

```jldoctest ex
julia> using Juliana;

julia> ct, structures, presc, fields, field_targets = Juliana.load_TROTS_data("Protons1");
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: slice.RescaleIntercept != -1000
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:142
┌ Warning: Not all CT slices have the exactly same spacing!
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/Dicom.jl:174
┌ Warning: WARNING: The parotid constraints aim to spare both, but the guideline demands only one.
└ @ Juliana /data/user/bellotti_r/differentiable-planning/src/TROTS.jl:148


```

Before we dive into the details of all the objects we just loaded, a word about units is needed.

32-bit floating point numbers are required unless stated otherwise.

## Units

All lengths in Juliana.jl are in units of ```cm``` and dose values are in ```Gy````, but this is a matter of scaling factors. Angles are expected to be in units of degrees.

## ScalarGrid

```julia
struct Grid{S<:Real, T<:Integer}
    spacing::SVector{3, S}
    origin::SVector{3, S}
    size::SVector{3, T}
end
```

```julia
struct ScalarGrid{A<:AbstractArray, S, T}
    data::A
    grid::Grid{S, T}
end
```

The CT image is of type ```Juliana.ScalarGrid```, which has a field ```data``` (a three-dimensional matrix)
and a field ```grid```, which is of type ```Juliana.Grid``` and contains the fields ```spacing, origin, size```. The fields of the grid are of type ```StaticArrays.SVector{3}```, i. e. their size is known at compile time.

The ```Juliana.ScalarGrid``` type is not limited to CTs. It can also be used to contain dose cubes, water-equivalent depth (WED) at the
nodes of a grid, binary masks etc. The latter is particularly useful to debug binary masks and distance calculations.

!!! warning Please notice that the HU-to-stopping-power lookup table assumes that the HU in the CT range from ```[-1000, 3096]```. A warning will be issued when loading a DICOM CT whose ``RescaleIntercept`` value is unequal to ``-1000``. This is simply to implement PSI conventions; feel free to ignore it for other lookup tables.

## Structure
```julia
struct Structure{A<:AbstractArray, B<:AbstractArray, C<:AbstractArray, S, T}
    name::String
    points::A
    mask::B
    distanceFromStructure::C
    grid::Grid{S, T}
    is_target::Bool
end
```

The field ```is_target``` indicates whether the structure is a target or not.

### Structure as a collection of points
A structure represents an anatomical region, such as targets, organs/at/risk (OARs) or other regions of interest.
It is identified by its ```name``` and the ```points``` of its contours. The latter is an array of shape ```n x 4```,
where ```n``` is the number of contour points. The first columns are the ```x, y, z``` coordinates of the contour in the
CT coordinate system.
The fourth column contains the contour index, which is semantically and integer but stored as a ```Float32```
due to the fact that all columns in arrays must have the same data type.

### Structure as a binary mask
Further, ```Juliana.Structure``` contains a so-called distance mask, which is an array of the same shape as the CT. Each voxel contains the closest distance of the voxel to any contour of the structure of interest. Positive distance means that the voxel lies outside the structure, negative distance that it lies inside the structure. A distance of zero means that the voxel lies on a contour border.

A binary mask is obtained via thresholding of the distance mask using ```Juliana.to_binary_mask()```. It is also an array of the same shape as the CT. A voxel value of ```0``` indicates that the voxel does not belong to the structure, a value of ```1``` indicates that it does belong to the structure.

!!! warning
    The field ```distanceFromStructure``` contains the in-slice distance by default. If there are no contour points in a slice, that slice's values are ```Inf```, which allows thresholding. If you want the full three-dimensional distance, you first have to calculate the in-slice distance and then call the function ```Juliana.calculate_outside_distance_3D()``` to take distance across z-slices also into account.

The ```grid``` property corresponds to the grid on which the ```mask``` and ```distanceFromStructure``` are defined.

## Utility functions for masks

### Mask to voxel coordinates

It is often needed to find the ```xyz``` coordinates for all voxels in a structure. For example, imagine you would like to calculate the 
maximum and minimum WED for all voxels within a structure. Such a task is useful for spot placement (energy selection) and involves calculating
the WED at each voxel (or a subset thereof) within the structure.

The function ```points = Juliana.mask_to_points(grid, mask)``` takes a three-dimensional array representing the binary mask of the structure and returns a matrix of size ```3 x n```, where the number ```n``` represents the number of voxels within the structure, i. e. ```n = sum(mask)```.
Each column in ```points``` represents the ```xyz``` coordinates of another voxel at which the mask is true. Like for the CT, voxel values represent values at the voxel centre.

### Mask to voxel coordinates and indices

If we would like to reconstruct a three-dimensional volume such as the full WED cube based on the voxel values at some points (see next section), for example for exporting the WED as a DICOM RT dose.
In such cases, not only the point's coordinates, but also the corresponding grid indices are needed.
This can be calculated using  ```points, indices = Juliana.mask_to_points_and_indices(grid, mask)```, where ```indices``` is also of shape ```3 x n```.
The ```xyz``` coordinats of the ```i```-th point are ```points[:, i]``` and its corresponding indices are ```indices[:, i]```.

## Points to cubes

A recurring task is to calculate some quantity of interest at a set of voxels and then fill this information into a cube for easy handling, visualisation and export to DICOM RT dose.
The utility function ```cube = Juliana.flat_vector_to_cube(grid, indices, values_flat)``` can do this. The resulting cube will have the values of the given flattened positions and zero everywhere else.

The short example in the next section walks you through the calculation of the WED for all voxels within the target while applying a checkerboard pattern and then transforms it back to a cube of the same size as the CT.

## Water-equivalent depth (wed)

We will now demonstrate how to calculate the water-equivalient depth (WED). First, we need to define the points at which to evaluate the WED.
This can be either a matrix of shape ```3 x n``` with each column representing an ```xyz``` coordinate at which to calculate the WED or a ```Juliana.Grid```. In the case of the matrix, the result type is a vector containing the WED values for each point. In the case of the matrix, the result will be a ```Juliana.ScalarGrid``` containing the WED at each grid node.

We will now demonstrate how to calculate the WED at the points within the target, while only considering the voxels based on a checkerboard pattern, and then converts it back to a ```Juliana.ScalarGrid```.


```jldoctest ex
target = structures["CTV High"]
checkerboard = Juliana.build_checker_board_mask(ct.grid)
mask = target.mask .&& checkerboard
points, indices = Juliana.mask_to_points_and_indices(ct.grid, mask);
```

The calculation of the WED requires us to load the ```Juliana.JulianaTps``` that contains the lookup curves for HU-to-stopping-power calculations.
We can calculate the densities to calculate the WED.

```jldoctest ex
tps = Juliana.load_juliana_tps()

gantry_angle = 34.f0 # degrees in [-30, 270] (?) for PSI
couch_angle = 0.f0 # degrees in [-180, 180] for PSI

densities = Juliana.ScalarGrid(
    tps.convert_to_sp.(ct.data),
    ct.grid
)

wed_at_points = Juliana.calculate_wed(
    densities,
    gantry_angle,
    couch_angle,
    points,
)

wed_cube = Juliana.flat_vector_to_cube(
    ct.grid,
    indices,
    wed_at_points,
);
```

Calculating the WED on the CT grid is even easier:

```jldoctest ex
wed_scalargrid = Juliana.calculate_wed(
    densities,
    gantry_angle,
    couch_angle,
    ct.grid,
);
```

WED calculation is parallelised on the GPU using ```KernelAbstractions.jl```. Currently, only CUDA is supported as the backend.

## Evaluating a cube at indices

```jldoctest ex
vector_flat = Juliana.volume_at_indices(
    ct,
    indices,
);
```

## In-beam mask
Often, you would like to identify the region covered by a field/beam. This can be calculated efficiently (parallelised on the GPU; only CUDA supported currently). The calculation of the binary mask is implemented as a ray casting algorithm.

```jldoctest ex
in_beam_mask = calculate_in_beam_mask(
    target.mask,
    gantry_angle,
    couch_angle,
    ct.grid,
);
```

## Angles to direction vector
```jldoctest ex
# direction is a vector of size 3. It points from the iso-centre towards the nozzle.
direction = Juliana.angles_to_direction(gantry, couch);
```

## Coordinate transformations

There are utility functions to go from ```xyz``` to ```stu``` coordinates and
versa. This can be done for single points or multiple points
(given as a matrix of size ```3 x n```) and is particularly useful for implementing
spot placement.

```jldoctest ex
using StatsBase
iso_center = mean(target.points, dims=1)[1:3]

points_stu = Juliana.XyzToStu(
    points,
    iso_center,
    gantry_angle,
    couch_angle,
)

points_xyz = Juliana.StuToXyz(
	points,
    iso_center,
    gantry_angle,
    couch_angle,
)
```

## Treatment plans
```julia
struct TreatmentPlan{T<:Real}
    fields::Vector{FieldDefinition{T}}
    prescribedDose::T
end


struct FieldDefinition{T<:Real}
    label::String
    id::Integer
    gantryAngle::T
    couchAngle::T
    nozzleExtraction::T
    fieldCenter::Dict{String, T} # keys: "x", "y", "z". Leftover from FIonA standalone.
    spots::Vector{Spot{T}}
end


struct Spot{T<:Real}
    id::Integer
    t::T
    u::T
    weight::T
    energykeV::T
    # 0 or 1
    numberOfAbsorbers::Int32
end
```

## Autoplanning

This simple function call optimises spot weights using the JulianA algorithm.

DOI: 10.48550/ARXIV.2305.10211

```jldoctest ex
plan_optimised = Juliana.optimise_juliana(
    tps,
    ct,
    structures,
    prescriptions,
    plan,
);
```

## Dose calculation

```jldoctest ex
resolution = 0.35f0 # cm

dose = Juliana.calculate_dose(
    tps,
    resolution,
    ct,
    plan,
);
```


## Extracting dose metrics (mean, max, DVH)

```jldoctest ex
target_dose = Juliana.hottest_target(presc)[2]

maximum_dose = maximum(target.mask * dose.data)

mean_dose = Juliana.mean_dose(dose, target.mask)

# Volume that receives at least 95% of the prescribed dose.
V95 = Juliana.dvh_v(dose.data, target.mask, 0.95 * target_dose)
# Dose receive by the coldest 10% of the target.
D90 = Juliana.dvh_d(dose.data, target.mask, 0.9f0);
```


## DICOM

```jldoctest ex
patient_ID = "my_name"
study_instance_uid = Juliana.get_study_instance_uid(patient_ID)
output_dir = "."

Juliana.dicom_export_to_directory(
    ct,
    structures,
    output_dir,
    study_instance_uid,
    patient_ID,
    Dict{String, Juliana.ScalarGrid}(
        "my_dose" => dose,
    ),
);
```

## Plotting

## Robust optimisation
Some code is there, but there has not been much time to test it, unfortunately.

## Spot position

# TPS interface

# FionaStandalone.jl


# Reference

```@autodocs
Modules = [Juliana]
```
