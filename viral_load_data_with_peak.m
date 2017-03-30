function [pat_num, time, viral_load, LOD] = viral_load_data_with_peak(serotype, type_infection, CM)

M = xlsread('Clapham_data.xlsx');

 
j = 1; 

for i = 1:size(M, 1)
    if M(i,2) == serotype && M(i,3) == type_infection && M(i,4) == CM  && M(i,1)~= M(i-1,1)
            pat_num_full(j) = M(i,1);
            j = j+1; 
    end
end

temp_viral_load = NaN(length(pat_num_full), 16);
time = NaN(length(pat_num_full), 16); 
viral_load = NaN(length(pat_num_full), 16);
LOD = NaN(length(pat_num_full), 16); 

k = 1; 


for i = 1:length(pat_num_full)
    temp = find(M(:,1) == pat_num_full(i)); 
    
     temp_viral_load(i,1:length(temp)) = log10(M(temp, 6));
     temp2 = find(temp_viral_load(i,:) == max(temp_viral_load(i,:)));
    
    if temp2 ~= 1 && isnan(temp_viral_load(i,1)) == 0 %&& (temp_viral_load(i, temp2(1)) - temp_viral_load(i,1)) > 0.2
        
        time(k,1:length(temp)) = M(temp, 5); 
        viral_load(k,1:length(temp)) = log10(M(temp, 6));
        LOD(k,1:length(temp)) = log10(M(temp, 7)); 
        pat_num(k) = pat_num_full(i); 

            for j = 1:length(temp)
                if viral_load(k,j) < LOD(k,j)
                    viral_load(k,j) = LOD(k,j); 
                end
            end
            
        k = k + 1;     
        
    end
end

time = time(1:length(pat_num),:); 
viral_load = viral_load(1:length(pat_num), :); 
LOD = LOD(1:length(pat_num), :); 
end