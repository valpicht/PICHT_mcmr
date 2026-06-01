"""
Defines the types and main methods for user-defined geometries.
"""
module Obstructions
include("fields.jl")
include("obstruction_types.jl")
include("obstruction_groups.jl")
import StaticArrays: MVector
import .ObstructionTypes: ObstructionType, fields
import .Fields: Field, FieldValue, isglobal, description, value_as_vector
import .ObstructionGroups: ObstructionGroup, IndexedObstruction, key_value_pairs, nvolumes
import ....Methods: get_rotation

function field_to_docs(key::Symbol, field_value::FieldValue{T}) where {T}
    key_txt = replace(String(key), "_"=>"\\_")
    return "- $(key_txt): $(field_to_docs(field_value))\n"
end

function field_to_docs(field_value::FieldValue{T}) where {T}
    field = field_value.field
    as_string = description(field_value)
    if ~isnothing(field.default_value)
        as_string = as_string * "  default value: $(field.default_value)"
    end
    return as_string
end

for obstruction_type in (
    ObstructionType(:Wall; ndim=1, volumes=[]),
    ObstructionType(
        :Cylinder; ndim=2, fields=[
            Field{Float64}(:radius, "Radius of the cylinder.", required=true), 
            Field{Float64}(:g_ratio, "Inner/outer radius used for susceptibility calculation", 1.),
            Field{Float64}(:susceptibility_iso, "Isotropic component of the susceptibility (in ppm um). It can be calculated by multiplying the isotropic susceptibility of the simulated tissue by the thickness.", -0.1),
            Field{Float64}(:susceptibility_aniso, "Anisotropic component of the susceptibility (in ppm um). It can be calculated by multiplying the anisotropic susceptibility of the simulated tissue by the thickness", -0.1),
            Field{Float64}(:lorentz_radius, "Only compute field explicitly for cylinders with this Lorentz radius.", 5.),
        ]),
    ObstructionType(
        :Annulus; plural=:Annuli, ndim=2, surfaces=[:inner_surface, :outer_surface], volumes=[:inner_volume, :outer_volume], fields=[
            Field{Float64}(:inner, "Radius of the inner cylinder.", required=true), 
            Field{Float64}(:outer, "Radius of the outer cylinder.", required=true), 
            Field{Bool}(:myelin, "Whether the annulus is myelinated.", false, required=true), 
            Field{Float64}(:susceptibility_iso, "Isotropic component of the myelin susceptibility (in ppm).", -0.1),
            Field{Float64}(:susceptibility_aniso, "Anisotropic component of the myelin susceptibility (in ppm).", -0.1),
            Field{Float64}(:lorentz_radius, "Only compute field explicitly for a annuli with this Lorentz radius.", 5.),
        ]),
    # Sphere geometry definition.
    # Modified to support SWC-derived substrates represented as overlapping spheres.
    ObstructionType(
        :Sphere; ndim=3, fields=[
            Field{Float64}(:radius, "Radius of the sphere.", required=true),

            # Allows multiple spheres to intersect, which is useful when reconstructing
            # a cellular substrate from SWC points and radii.
            Field{Bool}(:overlapping, "Whether overlapping spheres should be treated as permeable in the overlapping region.", false),
        ]),
    
    ObstructionType(
        :Triangle; plural=:Mesh, ndim=3, include_shift=false, group_volumes=true, fields=[
            Field{MVector{3, Int}}(:triangles, "Each triangle is defined by 3 vertices into the mesh.", required=true),
            Field{Vector{MVector{3, Float64}}}(:vertices, "Positions of the corners of the triangular mesh.", required=true, only_group=true),
            Field{Bool}(:myelin, "Whether the mesh is myelinated.", false, required=true), 
            Field{Float64}(:susceptibility_iso, "Isotropic component of the susceptibility (in ppm um). It can be calculated by multiplying the isotropic susceptibility of the simulated tissue by the thickness.", -0.1),
            Field{Float64}(:susceptibility_aniso, "Anisotropic component of the susceptibility (in ppm um). It can be calculated by multiplying the anisotropic susceptibility of the simulated tissue by the thickness.", -0.1),
            Field{Float64}(:lorentz_radius, "Only compute field explicitly for triangles with this Lorentz radius.", 5.),
            Field{Int}(:components, "Which component this triangle belongs to. If not provided explicitly, this will be determined based on connectivity.", required=false),
        ]),
    ObstructionType(
        :Ring; plural=:BendyCylinder, ndim=3, include_shift=false, fields=[
            Field{MVector{3, Float64}}(:control_point, "Control points defining the path of the cylinder.", required=true), 
            Field{Float64}(:radius, "Radius at each control point.", required=true), 
            Field{Int}(:nsamples, "Number of mesh vertices along each ring.", 100, only_group=true),
            Field{MVector{3, Int}}(:closed, "After how many repeats in each dimension does the cylinder connect with itself. If not set the cylinder is not closed.", only_group=true), 
            Field{Int}(:spline_order, "Sets the order of the b-spine interpolating between control points.", 3, required=true, only_group=true), 
            Field{Bool}(:myelin, "Whether the cylinder is myelinated.", false, required=true, only_group=true), 
            Field{Float64}(:susceptibility_iso, "Isotropic component of the susceptibility (in ppm um). It can be calculated by multiplying the isotropic susceptibility of the simulated tissue by the thickness.", -0.1),
            Field{Float64}(:susceptibility_aniso, "Anisotropic component of the susceptibility (in ppm um). It can be calculated by multiplying the anisotropic susceptibility of the simulated tissue by the thickness.", -0.1),
            Field{Float64}(:lorentz_radius, "Only compute field explicitly for triangles with this Lorentz radius.", 5.),
        ]),
)
    name_string = String(obstruction_type.singular)
    plural_name = String(obstruction_type.plural)
    (unique_keys, field_values) = key_value_pairs(obstruction_type, 0)
    unique_field_values = [field_values[key] for key in unique_keys]
    sort!(unique_field_values, by=fv->~fv.field.required)
    fields_string = join(field_to_docs.(unique_keys, unique_field_values))
    @eval const $(obstruction_type.singular) = IndexedObstruction{$(QuoteNode(obstruction_type.plural))}
    @eval const $(obstruction_type.plural) = ObstructionGroup{$(QuoteNode(obstruction_type.plural))}
    constructor = obstruction_type.plural
    single_type = obstruction_type.singular
    @eval begin
        $(constructor)(; kwargs...) = ObstructionGroup($(obstruction_type); kwargs...)
        @doc """
            $($constructor))(; fields...)

        Creates a set of [`MCMRSimulator.""" * $(name_string) * """`](@ref MCMRSimulator.Geometries.User.Obstructions.""" * $(name_string) * """) objects.
        Fields can be set using keyword arguments.
        The following fields are available:
        """ * $(fields_string) $constructor

        @doc """
        Singular """ * $(name_string) * """ object that is obtained by indexing a [`""" * $(plural_name) * """`](@ref) object.
        """ $single_type 
    end
end

end