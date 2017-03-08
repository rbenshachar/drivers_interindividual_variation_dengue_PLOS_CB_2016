function [pat_num, time, viral_load, LOD] = viral_load_data(serotype, type_infection, CM)

%Inputs: 
%serotype: 1 (DENV-1), 2 (DENV-2), 3 (DENV-3)
%type_infection: 1 (primary), 2 (secondary)
%CM: clinical manifestation: 1 (DF), 2 (DHF)

%Outputs:
%patnum:the patient IDs 
%three kinetics of the patient viral load measurements: 
%time of measurements (time), viral
%load measurements (viral_load) and limit of detection (LOD)

%modified_data comes from data from Clapham et al 2014, :Within-host viral dynamics of
%dengue serotype 1 infection," J R. Soc. Inteface
%modified diagnoses into ordinal variables for easy reading into MATLAB
%deleted columns for genotype, serology,sample_no, date_time,
%illness_time_at_admission_days, illness_time_at_admission_hrs,
%time_since_first_sample, time_next_sample, Id_code, qpcr_interpet
%serotype: 0 = neg, 1 = DENV-1, 2 = DENV-2, 3 = DENV-3, 4 = DENV-4
%primary_vs_secondary: 0 = neg, 1 = primary, 2= secondary
%diagnosis_WHO: 0 = neg, 1 = DF, 2 = DHF (any grade)
%time since illness onset - readin as time in days, not hours
%vl_ml_plasma and ld_vl_ml_plasma- rename NA, FUS and BA as blanks

M = xlsread('modified_data.xlsx');

j = 1; 

for i = 1:size(M, 1)
    if M(i,2) == serotype && M(i,3) == type_infection && M(i,4) == CM  && M(i,1)~= M(i-1,1) %making sure not counting same person twice
            pat_num(j) = M(i,1); %patient id
            j = j+1; 
    end
end

%initializing vectors for time, viral load and LOD for each patient
time = NaN(length(pat_num), 16); 
viral_load = NaN(length(pat_num), 16);
LOD = NaN(length(pat_num), 16); 


for i = 1:length(pat_num)
    temp = find(M(:,1) == pat_num(i)); 
    time(i,1:length(temp)) = M(temp, 5); 
    viral_load(i,1:length(temp)) = log10(M(temp, 6)); 
    LOD(i,1:length(temp)) = log10(M(temp, 7)); 
    
        %setting viral load to be LOD value if lower than LOD
        for j = 1:length(temp)
            if viral_load(i,j) < LOD(i,j)
                viral_load(i,j) = LOD(i,j); 
            end
        end
end

end