%save_params

%saving paramaters used in 'params' structure

params.time_start = 0; 
params.time_end = 35; 

%Set parameters - from TABLE 1 in PLOS CB paper
params.deltaT = 1e-6; 
params.alpha = 1e-3; 
params.omega = 1e4; 
params.d = 0.07; %parameter dN
params.dT = 0.1; 

%IC
params.Xinit = 1e7; 
params.Yinit = 0; 
params.Ninit = 0; 
params.Tinit = 1e5; 

%saving viral load data
[~, params.time_PI_1, params.data_PI_1, params.LOD_PI_1] = viral_load_data(1,1,1);
[~, params.time_PI_2, params.data_PI_2, params.LOD_PI_2] = viral_load_data(2,1,1);
[~, params.time_PI_3, params.data_PI_3, params.LOD_PI_3] = viral_load_data(3,1,1);

[~, params.time_SI_1_DF, params.data_SI_1_DF, params.LOD_SI_1_DF] = viral_load_data(1,2,1);
[~, params.time_SI_2_DF, params.data_SI_2_DF, params.LOD_SI_2_DF] = viral_load_data(2,2,1);
[~, params.time_SI_3_DF, params.data_SI_3_DF, params.LOD_SI_3_DF] = viral_load_data(3,2,1);

[~, params.time_SI_1_DHF, params.data_SI_1_DHF, params.LOD_SI_1_DHF] = viral_load_data(1,2,2);
[~, params.time_SI_2_DHF, params.data_SI_2_DHF, params.LOD_SI_2_DHF] = viral_load_data(2,2,2);
[~, params.time_SI_3_DHF, params.data_SI_3_DHF, params.LOD_SI_3_DHF] = viral_load_data(3,2,2);

%by CM 
params.time_PI = cat(1, params.time_PI_1, params.time_PI_2, params.time_PI_3); 
params.data_PI = cat(1, params.data_PI_1, params.data_PI_2, params.data_PI_3); 
params.LOD_PI = cat(1, params.LOD_PI_1, params.LOD_PI_2, params.LOD_PI_3); 

params.time_SI_DF = cat(1, params.time_SI_DF_1, params.time_SI_DF_2, params.time_SI_DF_3); 
params.data_SI_DF = cat(1, params.data_SI_DF_1, params.data_SI_DF_2, params.data_SI_DF_3); 
params.LOD_SI_DF = cat(1, params.LOD_SI_DF_1, params.LOD_SI_DF_2, params.LOD_SI_DF_3); 

params.time_SI_DHF = cat(1, params.time_SI_DHF_1, params.time_SI_DHF_2, params.time_SI_DHF_3); 
params.data_SI_DHF = cat(1, params.data_SI_DHF_1, params.data_SI_DHF_2, params.data_SI_DHF_3); 
params.LOD_SI_DHF = cat(1, params.LOD_SI_DHF_1, params.LOD_SI_DHF_2, params.LOD_SI_DHF_3); 

% by immune status
params.time_SI = cat(1, params.time_SI_DF, params.time_SI_DHF); 
params.data_SI = cat(1, params.data_SI_DF, params.data_SI_DHF); 
params.LOD_SI = cat(1, params.LOD_SI_DF, params.LOD_SI_DHF); 

save('params', 'params')