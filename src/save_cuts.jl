using JSON

function cut_to_dict(cut::DualSDDP.PrimalCut)
    return Dict(
        "value" => cut.value,
        "slope" => cut.slope,
        "x0" => cut.x0
    )
end

function cut_to_dict(cut::DualSDDP.DualCut)
    return Dict(
        "value" => cut.value,
        "slope_π" => cut.slope_π,
        "slope_γ" => cut.slope_γ,
        "π0" => cut.π0,
        "γ0" => cut.γ0
    )
end

function write_cuts_to_file(
    model,
    filename::String
) where T
    cuts = Dict{String, Any}[]
    for (i, stage) in enumerate(model)
        node_cuts = Dict(
            "node" => i,
            "single_cuts" => Dict{String, Any}[]
        )
        for cut in stage.ext[:cuts]
            push!(node_cuts["single_cuts"], cut_to_dict(cut))
        end
        push!(cuts, node_cuts)
    end

    open(filename, "w") do io
        JSON.print(io, cuts, 2)
    end
    return
end