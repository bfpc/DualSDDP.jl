import NPZ

eafs = NPZ.npzread("../eafs.npz")["eafs"]
eafs = eafs[:,:,:]
