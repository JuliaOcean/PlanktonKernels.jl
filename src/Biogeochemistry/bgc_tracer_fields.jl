const bgc_tracer_names = (:DIC, :NH4, :NO3, :PO4, :DFe, :O2, :DOC, :DON, :DOP,
                          :PIFe, :POC, :PON, :POP, :POFe, :Dust)

"""
    bgc_tracer_default_init()
Generate default bgc tracer initial conditions.
"""
function bgc_tracer_default_init()
    DIC =(init = 20.0, rand_noise = 0.1)
    NH4 =(init = 0.5, rand_noise = 0.1)
    NO3 =(init = 0.8, rand_noise = 0.1)
    PO4 =(init = 0.10, rand_noise = 0.1)
    DFe =(init = 1.0e-6, rand_noise = 0.1)
    O2 =(init = 160.0, rand_noise = 0.1)
    DOC =(init = 10.0, rand_noise = 0.1)
    DON =(init = 0.1, rand_noise = 0.1)
    DOP =(init = 0.05, rand_noise = 0.1)
    PIFe=(init = 0.0, rand_noise = 0.1)
    POC =(init = 0.0, rand_noise = 0.1)
    PON =(init = 0.0, rand_noise = 0.1)
    POP =(init = 0.0, rand_noise = 0.1)
    POFe=(init = 0.0, rand_noise = 0.1)
    Dust =(init = 0.0, rand_noise = 0.1)
    return (DIC=DIC, NH4=NH4, NO3=NO3, PO4=PO4, DFe=DFe, O2=O2, DOC=DOC, DON=DON,
            DOP=DOP, PIFe=PIFe, POC=POC, PON=PON, POP=POP, POFe=POFe,
            Dust=Dust)
end

"""
    generate_bgc_tracers(arch, grid, source=bgc_tracer_default_init(), FT=Float32)
Initialize the standard tracers through Fields, preserving their canonical order.
Source maps tracer names to `(init=value, rand_noise=noise)` or serialized array
paths, as accepted by Fields.generate_tracers. All standard names are required.
"""
function generate_bgc_tracers(arch::Architecture, grid::AbstractGrid,
                          source::NamedTuple=bgc_tracer_default_init(), FT::DataType=Float32)
    return Fields.generate_tracers(arch, grid, bgc_tracer_names, source, FT)
end
