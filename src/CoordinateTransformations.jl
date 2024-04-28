# using StaticArrays


# function xyz_to_stu_matrix(gantry_angle, couch_angle)
#     θ = gantry_angle
#     φ = -couch_angle

#     # The (1, -1, -1) is only for head-first-supine position, which is used
#     # at PSI. (?)
#     M = @SMatrix [
#           sin(θ)*cos(φ)  sin(θ)*sin(φ)   cos(θ)
#           cos(θ)*cos(φ)  cos(θ)*sin(φ)  -sin(θ)
#         -sin(θ)               cos(φ)       0
#     ]
#     return M .* (1, -1, -1)
# end


"""
    XyzToStu(px, py, pz, XyzIso, gantry_angle::Float32, couch_angle::Float32)   

This function does the conversion from the xyz coordinates (of the CT) to the
stu coordinates (the ones of the beam). gantry_angle and couch_angle must be
in degrees.
"""
function XyzToStu(px, py, pz, XyzIso, gantry_angle::Float32, couch_angle::Float32)
    Delta1 = px - XyzIso[1]
    Delta2 = py - XyzIso[2]
    Delta3 = pz - XyzIso[3]

    # We do not use the xyz_to_stu matrix because that matrix-vector multiplication
    # allocates memory, which is not possible on the GPU.
    s = - sind(gantry_angle)*cosd(couch_angle)*Delta1 + cosd(gantry_angle)*Delta2 - sind(gantry_angle)*sind(couch_angle)*Delta3
    t =   cosd(gantry_angle)*cosd(couch_angle)*Delta1 + sind(gantry_angle)*Delta2 + cosd(gantry_angle)*sind(couch_angle)*Delta3
    u =                      sind(couch_angle)*Delta1                             -                    cosd(couch_angle)*Delta3

    return s, t, u
end


function XyzToStu(px::Vector{T}, py::Vector{T}, pz::Vector{T}, XyzIso, gantry_angle::Float32, couch_angle::Float32) where {T <: Real}
    stu = Matrix{T}(undef, 3, length(px))
    for i in 1:length(px)
        stu[:, i] .= XyzToStu(
            px[i],
            py[i],
            pz[i],
            XyzIso,
            gantry_angle,
            couch_angle,
        )
    end
    return stu
end


function XyzToStu(points::Matrix{T}, XyzIso, gantry_angle::Float32, couch_angle::Float32) where {T <: Real}
    return XyzToStu(
        points[1, :],
        points[2, :],
        points[3, :],
        XyzIso,
        gantry_angle,
        couch_angle,
    )
end


function XyzToStu(points, field::Juliana.FieldDefinition)
    return Juliana.XyzToStu(
        points,
        (
            field.fieldCenter["x"],
            field.fieldCenter["y"],
            field.fieldCenter["z"],
        ),
        field.gantryAngle,
        field.couchAngle,
    )
end


"""
    StuToXyz(s, t, u, XyzIso, gantry_angle::Float32, couch_angle::Float32)   

This function does the conversion from the stu coordinates to the
xyz coordinates. gantry_angle and couch_angle must be in degrees.

"""
function StuToXyz(s, t, u, XyzIso, gantry_angle::Float32, couch_angle::Float32)
    # The inverse of a rotation is the rotation with negative angles.
    # + xyzIso: Needed because the implementation transforms the delta, Not
    # the point. (see XyzToStu)
    # TODO: Solve this in a cleaner way.
    x = XyzIso[1] - sind(gantry_angle)*cosd(couch_angle)*s + cosd(gantry_angle)*cosd(couch_angle)*t + sind(couch_angle)*u
    y = XyzIso[2] + cosd(gantry_angle)                  *s + sind(gantry_angle)                  *t
    z = XyzIso[3] - sind(gantry_angle)*sind(couch_angle)*s + cosd(gantry_angle)*sind(couch_angle)*t - cosd(couch_angle)*u

    return x, y, z
end


function StuToXyz(s::Vector{T}, t::Vector{T}, u::Vector{T}, XyzIso, gantry_angle::Float32, couch_angle::Float32) where {T <: Real}
    xyz = Matrix{T}(undef, 3, length(s))
    for i in 1:length(s)
        xyz[:, i] .= StuToXyz(
            s[i],
            t[i],
            u[i],
            XyzIso,
            gantry_angle,
            couch_angle,
        )
    end
    return xyz
end


function Juliana.StuToXyz(points::Matrix{T}, XyzIso, gantry_angle::Float32, couch_angle::Float32) where {T <: Real}
    return Juliana.StuToXyz(
        points[1, :],
        points[2, :],
        points[3, :],
        XyzIso,
        gantry_angle,
        couch_angle,
    )
end


function Juliana.StuToXyz(points::Matrix{T}, field::Juliana.FieldDefinition) where {T <: Real}
    iso_center = (
        field.fieldCenter["x"],
        field.fieldCenter["y"],
        field.fieldCenter["z"],
    )
    return Juliana.StuToXyz(
        points[1, :],
        points[2, :],
        points[3, :],
        iso_center,
        field.gantryAngle,
        field.couchAngle,
    )
end
