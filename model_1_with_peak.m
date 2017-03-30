%model_1_with_peak

%model 1 fit to peak viral load dataset

clear all; 
close all; 
clc; 

load('params');

iters = 300000; 
l = 200000; 
m = 100;   

data_with_peak;

startvalue = set_IC(1,[1e-10, 5, 1e-3, 5, 0.1, 1e-6]);
sd = [1e-10,  0.1, 1e-4, 0.1, 0.01, 1e-6];

[chain, accep_prop, my_posterior] = run_metropolis_full_peak(startvalue, iters, params,1, sd); 
 
save('chain_1_IC_1_peak', 'chain');
save('my_posterior_1_IC_1_peak', 'my_posterior');
save('accep_prop_1_IC_1_peak', 'accep_prop') 
 
