import NPZ

fname = joinpath(@__DIR__, "..", "eafs.npz")
eafs = NPZ.npzread(fname)["eafs"]
eafs = eafs[:,:,:]
