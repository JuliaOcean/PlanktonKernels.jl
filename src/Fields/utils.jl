##### adapt the field related structs to GPU
Adapt.adapt_structure(to, field::Field{FT}) where {FT} =
    Field{FT}(
        Adapt.adapt(to, field.data),
        Adapt.adapt(to, field.bc))

Adapt.adapt_structure(to, bc::BoundaryCondition) =
    BoundaryCondition(
        Adapt.adapt(to, bc.c),
        Adapt.adapt(to, bc.v))

Adapt.adapt_structure(to, bcs::BoundaryConditions) =
    BoundaryConditions(
        Adapt.adapt(to, bcs.west),
        Adapt.adapt(to, bcs.east),
        Adapt.adapt(to, bcs.north),
        Adapt.adapt(to, bcs.south),
        Adapt.adapt(to, bcs.top),
        Adapt.adapt(to, bcs.bottom))