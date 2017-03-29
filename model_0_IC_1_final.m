%model_0_IC_1_final

%run MCMC algorithm for model 0 from PLOS CB 2016 paper
%model 0

clear all; 
close all; 
clc; 

load('params');

iters = 300000; 
l = 200000; 
m = 100;   

sd = [1e-10,  0.1, 1e-4,0.01, 0.1];

startvalue = set_IC(1,[1e-10, 5,1e-4, 0.1,-2]); 

[chain, accep_prop, my_posterior] = run_metropolis_full(startvalue, iters, params,0, sd, 26, 202); 
 
save('chain_0_IC_1_final', 'chain');
save('my_posterior_0_IC_1_final', 'my_posterior');
save('accep_prop_0_IC_1_final', 'accep_prop') 
 
