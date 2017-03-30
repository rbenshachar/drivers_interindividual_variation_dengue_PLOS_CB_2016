%OAS_IC_1_final

%OAS model

%difference in deltaT by secondary infection clinical manifestation

clear all; 
close all; 
clc; 

load('params');

iters = 300000; 
l = 200000; 
m = 100; 

startvalue = set_IC(1, [5e-10, 5, 1e-4, 0.1, 1e-6, 1e-5, -2]);

sd = [1e-10,  0.1, 1e-4, 0.01, 1e-7, 1e-9, 0.1]; 

[chain_CM, accep_prop_CM, my_posterior_CM] = run_metropolis_MCMC_CM(startvalue, iters, params,1, sd); 

save('chain_OAS_IC_1_final', 'chain_CM');
save('my_posterior_OAS_IC_1_final', 'my_posterior_CM');
save('accep_prop_OAS_IC_1_final', 'accep_prop_CM') 

 
