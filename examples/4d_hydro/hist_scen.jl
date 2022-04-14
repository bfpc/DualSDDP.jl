import NPZ

fname = joinpath(@__DIR__, "..", "eafs.npz")
eafs = NPZ.npzread(fname)["eafs"]
if Main.nscen > 0
  maxscen = size(eafs, 2)
  k = min(maxscen, Main.nscen)
  eafs = eafs[:,1:k,:]
end
