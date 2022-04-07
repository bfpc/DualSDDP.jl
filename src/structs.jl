# All are function of time x branch 
struct MSLBO
  A :: Function # Ax_t + Bx_{t-1} +Ty = d
  B :: Function 
  T :: Function
  c :: Function # marginal cost of y_t
  d :: Function
  Ux :: Function #upper bound on the state x
  Uy :: Function #upper bound on the control y
  lb :: Function # lower bound on the value of the problem
  ub :: Function # upper bound on the value of the problem
  Lip :: Function # upper bound on the Lipschitz constant
  prob :: Function # reference probability over branches
end
