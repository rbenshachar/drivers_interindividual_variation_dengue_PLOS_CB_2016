%ADE_IC_1_with_peak

%ADE model fit to subset of data with peak

%CM difference in beta 

clear all; 
close all; 
clc; 

load('params');

iters = 300000; 
l = 200000; 
m = 100; %100;  
 
data_with_peak;

sd = [1e-10,  0.1, 1e-4, 0.1, 0.01, 1e-7, 1e-12]; 

startvalue = set_IC(1, [5e-10, 5, 1e-4, 5, 0.1, 1e-6, 1e-11]);
 
[chain_CM, accep_prop_CM, my_posterior_CM] = run_metropolis_MCMC_CM_peak(startvalue, iters, params,2, sd); 

save('chain_ADE_IC_1_peak', 'chain_CM');
save('my_posterior_ADE_IC_1_peak', 'my_posterior_CM');
save('accep_prop_ADE_IC_1_peak', 'accep_prop_CM') 
 
