import NPZ

eafs = NPZ.npzread("../eafs.npz")["eafs"]
eafs = eafs[:,1:40,:]
