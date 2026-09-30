const bgc_tracer_names = (:DIC, :NH4, :NO3, :PO4, :DFe, :O2, :DOC, :DON, :DOP,
                      :PFe_inorg, :POC, :PON, :POP, :PFe_bio, :Dust)

"""
    bgc_tracer_init()
Generate defalut bgc tracer initial conditions.
"""
function bgc_tracer_init()
    init = (DIC=20.0, NH4=0.5, NO3=0.8, PO4=0.10, DFe=1.0e-6, O2=160.0, DOC=10.0, DON=0.1, DOP=0.05, PFe_inorg=0.0, POC=0.0, PON=0.0, POP=0.0, PFe_bio=0.0, Dust=0.0)
    rand_noise = (DIC=0.1, NH4=0.1, NO3=0.1, PO4=0.1, DFe=0.1, O2=0.1, DOC=0.1, DON=0.1, DOP=0.1, PFe_inorg=0.1, POC=0.1, PON=0.1, POP=0.1, PFe_bio=0.1, Dust=0.1)
    return (initial_condition = init, rand_noise = rand_noise)
end

"""
    generate_tracers(arch, grid, source=bgc_tracer_init(), FT=Float32)
Initialize the standard tracers through Fields, preserving their canonical order.
Source is an initial_condition/rand_noise NamedTuple or a dictionary of serialized
array paths, as accepted by Fields.generate_tracers. All standard names are required.
"""
function generate_bgc_tracers(arch::Architecture, grid::AbstractGrid,
                          source::Union{Dict,NamedTuple}=bgc_tracer_init(), FT::DataType=Float32)
    return Fields.generate_tracers(arch, grid, bgc_tracer_names, source, FT)
end
