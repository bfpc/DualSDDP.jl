function readblock(file; nheader::Int=2)
  for i = 1:nheader
    line = readline(file)
  end
  line = readline(file)
  nblocks = length(split(line))

  idxs = Int[]
  data = Vector{Float64}[]
  while true
    line = readline(file)
    if line == nothing || line == ""
      println("Unclean break: no termination line")
      return idxs, data
    end

    parts = split(line)
    # 999... termination
    if length(parts) == 1
      if all([x == '9' for x in parts[1]])
        line = readline(file)
        if line == nothing || line == ""
          return idxs, data
        else
          println("Unclean break: no blank line after termination")
          return idxs, data
        end
      else
        println("Strange line: only index data, not termination")
        continue
      end
    end

    n = parse(Int, parts[1])
    push!(idxs, n)
    push!(data, [parse(Float64, x) for x in parts[2:end]])
  end
end
