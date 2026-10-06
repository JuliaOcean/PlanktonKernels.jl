##### calculate simple remineralization of DOM and POM as well as a simple nitrification
function bgc_tracer_forcing!(F, bgc_tracer_tmp, tracer, params, ΔT)
    for name in bgc_tracer_names
        @inbounds bgc_tracer_tmp[name].data .= max.(0.0f0, tracer[name].data)
    end

    @inbounds F.DIC.data .= F.DIC.data .+ bgc_tracer_tmp.DOC.data .* params["kDOC"] .* ΔT

    @inbounds F.DOC.data .= F.DOC.data .- bgc_tracer_tmp.DOC.data .* params["kDOC"] .* ΔT .+
                                          bgc_tracer_tmp.POC.data .* params["kPOC"] .* ΔT

    @inbounds F.POC.data .= F.POC.data .- bgc_tracer_tmp.POC.data .* params["kPOC"] .* ΔT

    @inbounds F.NH4.data .= F.NH4.data .+ bgc_tracer_tmp.DON.data .* params["kDON"] .* ΔT .-
                                          bgc_tracer_tmp.NH4.data .* params["Nit"]  .* ΔT

    @inbounds F.NO3.data .= F.NO3.data .+ bgc_tracer_tmp.NH4.data .* params["Nit"]  .* ΔT

    @inbounds F.DON.data .= F.DON.data .- bgc_tracer_tmp.DON.data .* params["kDON"] .* ΔT .+
                                          bgc_tracer_tmp.PON.data .* params["kPON"] .* ΔT

    @inbounds F.PON.data .= F.PON.data .- bgc_tracer_tmp.PON.data .* params["kPON"] .* ΔT

    @inbounds F.PO4.data .= F.PO4.data .+ bgc_tracer_tmp.DOP.data .* params["kDOP"] .* ΔT

    @inbounds F.DOP.data .= F.DOP.data .- bgc_tracer_tmp.DOP.data .* params["kDOP"] .* ΔT .+
                                          bgc_tracer_tmp.POP.data .* params["kPOP"] .* ΔT

    @inbounds F.POP.data .= F.POP.data .- bgc_tracer_tmp.POP.data .* params["kPOP"] .* ΔT

    ##### Iron cycle

    @inbounds F.DFe.data .= F.DFe.data .+ bgc_tracer_tmp.POFe.data .* params["kPOC"] .* ΔT .-
                            (params["lambda_POC"] .* bgc_tracer_tmp.POC.data .+ params["lambda_min"] .+
                             params["lambda_dust"] .* bgc_tracer_tmp.Dust.data) .* bgc_tracer_tmp.DFe.data .* params["DFeFrac"] .* ΔT .+ 
                             bgc_tracer_tmp.PIFe.data .* params["kdiss"] .* ΔT .- 
                             max.((bgc_tracer_tmp.DFe.data .- params["ligand"]), 0.0f0) .* bgc_tracer_tmp.DFe.data .* params["DFeFrac"] .* params["lambda_Fe"] .* ΔT

    @inbounds F.POFe.data .= F.POFe.data .+ bgc_tracer_tmp.POC.data .* 
                             bgc_tracer_tmp.DFe.data .* params["DFeFrac"] .* params["lambda_POC"] .* ΔT .- 
                                params["kPOC"] .* bgc_tracer_tmp.POFe.data  .* ΔT

    @inbounds F.PIFe.data .= F.PIFe.data .- bgc_tracer_tmp.PIFe.data .*
                             params["kdiss"].* ΔT .+ (params["lambda_min"] .+ params["lambda_dust"] .* 
                             bgc_tracer_tmp.Dust.data) .* bgc_tracer_tmp.DFe.data .* params["DFeFrac"] .* ΔT .+
                             max.((bgc_tracer_tmp.DFe.data .- params["ligand"]), 0.0f0) .* bgc_tracer_tmp.DFe.data .* params["DFeFrac"] .* params["lambda_Fe"] .* ΔT


    return nothing
end

