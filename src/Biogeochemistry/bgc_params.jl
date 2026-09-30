"Default reaction and tracer-transport parameters from PlanktonIndividuals (rates per second)."
function bgc_params_default(FT=Float32)
    params = Dict{String, FT}(
        "kDOC"         => 1/30/86400,         # Remineralization rate for DOC, turn over time: a month (per second)
        "Nit"          => 1/30/86400,         # Nitrification rate for NH4
        "kDON"         => 1/30/86400,         # Remineralization rate for DON, turn over time: a month (per second)
        "kDOP"         => 1/30/86400,         # Remineralization rate for DOP, turn over time: a month (per second)
        "kPOC"         => 1/30/86400,         # Remineralization rate for POC, turn over time: a month (per second)
        "kPON"         => 1/30/86400,         # Remineralization rate for PON, turn over time: a month (per second)
        "kPOP"         => 1/30/86400,         # Remineralization rate for POP, turn over time: a month (per second)
        "κh"           => 0.0e-6,             # Horizontal diffusion
        "κv"           => 0.0e-6,             # Vertical diffusion
        "w_sink_org"   => 2.31e-8,            # Settling velocity of organic matter (m/s)
        "w_sink_inorg" => 2.31e-5,            # Settling velocity of inorganic matter (m/s)
        "lambda_POC"   => 5.8e-8,             # Scavenging rate of iron by POC (m³/mmol/second)
        "lambda_min"   => 3.5e-10,            # Minimum scavenging rate of iron (per second)
        "kdiss"        => 4.4e-7 ,            # Dissolution rate of inorganic particulate iron  (per second)
        "lambda_dust"  => 1.7e-3,             # Scavenging rate of iron by dust (m³/kg/second)
        "lambda_Fe"    => 1.7,                # Coagulation rate of dissolved iron (m³/mmol/second)
        "ligand"       => 6.0e-4,             # Ligand concentration (mmol/m³/second)
        "DFeFrac"      => 0.01,               # Free DFe Fraction
        )
    return params
end

"""
    update_bgc_params(tmp::Ditc, FT::DataType)
Update parameter values based on a `Dict` provided by user

Keyword Arguments
=================
- `tmp`: a `Dict` containing the parameters needed to be upadated
- `FT`: Floating point data type. Default: `Float32`.
"""
function update_bgc_params(tmp::Dict, FT::DataType)
    parameters = bgc_params_default(FT)
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: bgc parameter not found $key"))
        else
            parameters[key] = FT(tmp[key])
        end
    end
    return parameters
end