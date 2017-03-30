%SS_3_IC_1_final

%serotype-specific differences in qT

clear all; 
close all; 
clc; 

load('params');

iters = 300000; 
l = 200000; 
m = 100; 
 
startvalue = set_IC(1, [5e-10, 5, 1e-4, 0.1, 1e-6, 1e-6, 1e-6, -2]);

sd = [1e-10,  0.1, 1e-4, 0.01, 1e-7, 1e-7, 1e-7, 0.1]; 

 
[chain_multi_level, accep_prop_multi_level, my_posterior_multi_level] = run_multi_level_metropolis_MCMC(startvalue, iters, params,3, sd); 

save('chain_SS_3_IC_1_final', 'chain_multi_level');
save('my_posterior_SS_3_IC_1_final', 'my_posterior_multi_level');
save('accep_prop_SS_3_IC_1_final', 'accep_prop_multi_level')
 
