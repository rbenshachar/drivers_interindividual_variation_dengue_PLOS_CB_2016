%model_1_IC_1_final

%model 1

clear all; 
close all; 
clc; 

load('params');

iters = 300000; 
l = 200000; 
m = 100;   

%beta, kappa, q, sigma, qT, V0
sd = [1e-10,  0.1, 1e-4, 0.01, 1e-7, 0.1];

startvalue = set_IC(1,[1e-10, 5, 1e-4, 0.1, 1e-6, -2]);

[chain, accep_prop, my_posterior] = run_metropolis_full(startvalue, iters, params,1, sd, 26, 202); 
 
save('chain_1_IC_1_final', 'chain');
save('my_posterior_1_IC_1_final', 'my_posterior');
save('accep_prop_1_IC_1_final', 'accep_prop') 
 
