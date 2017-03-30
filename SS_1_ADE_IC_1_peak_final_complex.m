%SS_1_ADE_IC_1_peak_final_complex

%differences in beta by CM and by serotype fit to viral peak dataset

clear all; 
close all; 
clc; 

load('params');

iters = 300000; 
l = 200000; 
m = 100; 
 
data_with_peak;

startvalue = set_IC(1, [5e-10, 5, 1e-4, 5, 0.1, 1e-6, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11]);

sd = [1e-10,  0.1, 1e-4, 0.1, 0.01, 1e-7, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12]; 
 
[chain_multi_level, accep_prop_multi_level, my_posterior_multi_level] = run_multi_level_metropolis_MCMC_peak(startvalue, iters, params,5, sd); 

save('chain_SS_1_ADE_IC_1_peak_final_prior_complex', 'chain_multi_level');
save('my_posterior_SS_1_ADE_IC_1_peak_final_prior_complex', 'my_posterior_multi_level');
save('accep_prop_SS_1_ADE_IC_1_peak_final_prior_complex', 'accep_prop_multi_level')
 
