%data_with_peak

%identifying subset of individuals with viral peak and puts this data in params structure

[~, params.p_time_PI_1, params.p_data_PI_1, params.p_LOD_PI_1] = viral_load_data_with_peak(1,1,1);
[~, params.p_time_PI_2, params.p_data_PI_2, params.p_LOD_PI_2] = viral_load_data_with_peak(2,1,1);
[~, tparams.p_time_PI_3, params.p_data_PI_3, params.p_LOD_PI_3] = viral_load_data_with_peak(3,1,1); 

[~, params.p_time_SI_DF_1,params.p_data_SI_DF_1, params.p_LOD_SI_DF_1] = viral_load_data_with_peak(1,2,1);
[~, params.p_time_SI_DF_2, params.p_data_SI_DF_2, params.p_LOD_SI_DF_2] = viral_load_data_with_peak(2,2,1);
[~, params.p_time_SI_DF_3, params.p_data_SI_DF_3, params.p_LOD_SI_DF_3] = viral_load_data_with_peak(3,2,1);

[~, params.p_time_SI_DHF_1, params.p_data_SI_DHF_1, params.p_LOD_SI_DHF_1] = viral_load_data_with_peak(1,2,2);
[~, params.p_time_SI_DHF_1, params.p_data_SI_DHF_2, params.p_LOD_SI_DHF_2] = viral_load_data_with_peak(2,2,2);
[~, params.p_time_SI_DHF_1, params.p_data_SI_DHF_3, params.p_LOD_SI_DHF_3] = viral_load_data_with_peak(3,2,2);

params.p_time_SI_1 = cat(1, params.p_time_SI_DF_1, params.p_time_SI_DHF_1); 
params.p_data_SI_1 = cat(1, params.p_data_SI_DF_1, params.p_data_SI_DHF_1); 
params.p_LOD_SI_1 = cat(1, params.p_LOD_SI_DF_1, params.p_LOD_SI_DHF_1); 

params.p_time_SI_2 = cat(1, params.p_time_SI_DF_2, params.p_time_SI_DHF_2); 
params.p_data_SI_2 = cat(1, params.p_data_SI_DF_2, params.p_data_SI_DHF_2); 
params.p_LOD_SI_2 = cat(1, params.p_LOD_SI_DF_2, params.p_LOD_SI_DHF_2); 

params.p_time_SI_3 = cat(1, params.p_time_SI_DF_3, params.p_time_SI_DHF_3); 
params.p_data_SI_3 = cat(1, params.p_data_SI_DF_3, params.p_data_SI_DHF_3); 
params.p_LOD_SI_3 = cat(1, params.p_LOD_SI_DF_3, params.p_LOD_SI_DHF_3); 

params.p_time_PI = cat(1, time_PI_DF_1, time_PI_DF_2, time_PI_DF_3);
params.p_data_PI = cat(1, data_PI_DF_1, data_PI_DF_2, data_PI_DF_3);
params.p_LOD_PI = cat(1, LOD_PI_DF_1, LOD_PI_DF_2, LOD_PI_DF_3);  

params.p_time_SI_DF = cat(1, time_SI_DF_1, time_SI_DF_2, time_SI_DF_3);
params.p_data_SI_DF = cat(1, data_SI_DF_1, data_SI_DF_2, data_SI_DF_3);  
params.p_LOD_SI_DF = cat(1, LOD_SI_DF_1, LOD_SI_DF_2, LOD_SI_DF_3); 

params.p_time_SI_DHF = cat(1, time_SI_DHF_1, time_SI_DHF_2, time_SI_DHF_3); %, time_SI_DHF_4); 
params.p_data_SI_DHF = cat(1,data_SI_DHF_1, data_SI_DHF_2, data_SI_DHF_3); %, data_SI_DHF_4);
params.p_LOD_SI_DHF = cat(1, LOD_SI_DHF_1, LOD_SI_DHF_2, LOD_SI_DHF_3); %, LOD_SI_DHF_4)

params.p_time_SI = cat(1, params.p_time_SI_DF, params.p_time_SI_DHF); 
params.p_data_SI = cat(1, params.p_data_SI_DF, params.p_data_SI_DHF); 
params.p_LOD_SI = cat(1, params.p_LOD_SI_DF, params.p_LOD_SI_DHF); 

