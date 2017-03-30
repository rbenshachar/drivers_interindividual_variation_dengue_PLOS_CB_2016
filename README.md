# drivers_interindividual_variation_dengue_PLOS_CB_2016
MCMC code for 2016 PLOS Computational Biology paper, "Drivers of inter-individual variation in dengue viral load dynamics" by R Ben-Shachar, S Schmidler and K Koelle. 
All code is written in MATLAB

run save_params.m to get set parameters for 'params' structure                                 

FILES TO RUN EACH MODEL IN THE PAPER:                    

Models fit to full dataset:                           
Model 0:   run model_0_IC_1_final.m                           
Model 1:   run model_1_IC_1_final.m                 
Model OAS: run OAS_1_IC_1_final.m
Model ADE: run ADE_1_IC_1_final.m                            
Model SSbeta: run SS_1_IC_1_final.m
Model SSq: run SS_2_IC_1_final.m                  
Model SSqT: run SS_3_IC_1_final.m                      
Model SSbetaADE: run SS_1_ADE_IC_1_final_complex.m             

Models fit to peak viral load dataset:                            
Model 1: run model_1_with_peak.m                     
Model ADE: run ADE_IC_1_with_peak.m            
Model SSbeta:
Model SSq:
Model SSqT: 
Model SSbetaADE: 

FINAL MCMC CHAINS AND POSTERIORS FOR EACH MODEL: 
Model 0: chain_0_IC_1_final.mat
Model 1: chain_1_IC_1_final.mat, chain_1_IC_2_final.mat, chain_1_IC_3_final.mat, chain_1_IC_4_final.mat (four different initial conditions)          
Model OAS:
Model ADE:
Model SSbeta:
Model SSq:
Model SSqT:

SUPPLEMENTARY MODELS: 
Model 1 fit separately to subset of individuals who received chloroquine (CQ) and those that received placebo (P): 
chain_1_IC_1_final_CQ.mat, chain_1_IC_1_final_P.mat
 
Model 1 fit for different deltaT values: 
In model 1 described above, deltaT = 10^-6. Fits of model 1 with values of deltaT = 10^-4, 10^-5, 10^-7 and 10^-8 were also fit. chain_1_IC_1_deltaT_i.mat corresponds to model 1 run with deltaT = 10^-i, for i = 4, 5, 7, 8

TO REPRODUCE FIGURES IN THE PAPER: 
Figure 1: run plot_serotype_viral_load_data.m            
Figure 2: run plot_model_1_dynamics_with_observation_noise.m              
Figure 3: run plot_model_1_dynamics_with_variance.m              
Figure 4: run plot_CM_densities.m                  
Figure 5: run plot_SS_densities.m                   
Figure 6: run plot_best_SS_models.m (also creates Figure S6)                

Figure S1: run plot_model_traces.m                        
Figure S2: run plot_model_1_correlations.m                             
Figure S3: run plot_different_deltaT.m                        
Figure S4: run plot_model_1_CQ.m                  
Figure S5: run plot_ADE_dynamics_peak.m                     
Figure S6: (see plot_best_SS_models.m)                
Figure S7: run plot_SS_complex_densities.m             

TO REPRODUCE TABLES IN THE PAPER: (to be added)
