library(ICIKendallTau)
set.seed(1234)
x = rnorm(50000)
y = rnorm(50000)
ici_kt(x, y, output = "other")

# Expected output when it's working:
# min_x: -4.226628
# min_y: -4.608599
# n_entry: 50000
# tot: 1249975000
# sum_obs: 50001
# dis: 625759472
# con_minus_dis (k_numerator): -1543944.000000
# n_tie: 0.000000
# m: 2499950000
# x_tie: 0.000000
# y_tie: 0.000000
# s_adjusted: -1543944.000000
# var: 13889305541666.666016
# z_b: -0.414277
# tau: -0.001235
# tau_max:1.000000
# tau      pvalue     tau_max 
# -0.00123518  0.67867094  1.00000000 

# what we get when it's not:
# min_x: -4.226628
# min_y: -4.608599
# n_entry: 50000
# tot: -897508648 # this is the problem right here
# sum_obs: 50001
# dis: 625759472
# con_minus_dis (k_numerator): -2149027592.000000
# n_tie: 0.000000
# m: -1795017296
# x_tie: 0.000000
# y_tie: 0.000000
# s_adjusted: -897508648.000000
# var: -9972816927026.666016
# z_b: nan
# tau: -1.000000
# tau_max:-1.000000
# tau  pvalue tau_max 
# -1     NaN      -1 

set.seed(1234)
x = rnorm(10000)
y = rnorm(10000)
ici_kt(x, y, output = "other")
