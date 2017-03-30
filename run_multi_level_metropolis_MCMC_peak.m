 function [chain, accep_prop, my_posterior] = run_multi_level_metropolis_MCMC_peak(startvalue, iterations, params, k, sd)

 if k == 1 || k == 2 || k == 3
     temp_chain = zeros(iterations, 8); 
     temp_chain(1,:) = startvalue;
 elseif k == 4
     temp_chain = zeros(iterations, 9); 
     temp_chain(1,:) = startvalue;
 elseif k == 5
     temp_chain = zeros(iterations,11); 
     temp_chain(1,:) = startvalue;
 end
  
 time_PI_1 = params.p_time_PI_1; 
 data_PI_1 = params.p_data_PI_1; 
 LOD_PI_1 = params.p_LOD_PI_1; 
 
       %take out all Nans
    for i = 1:9
        index = isnan(data_PI_1(i,:));
        time_temp = time_PI_1(i,:);
        data_temp  = data_PI_1(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_PI_1(i).time = time_temp; 
        indiv_PI_1(i).data = data_temp; 
        indiv_PI_1(i).LOD =  LOD_PI_1(i,1);
        indiv_PI_1(i).c = zeros(1, length(indiv_PI_1(i).data)); 
             for j = 1:length(indiv_PI_1(i).data)
                        if data_temp(j) <= indiv_PI_1(i).LOD      
                             indiv_PI_1(i).c(j) = 0;  
                        else
                             indiv_PI_1(i).c(j)  = 1; 
                        end
            end
    end
 
 time_PI_2 = params.p_time_PI_2; 
 data_PI_2 = params.p_data_PI_2; 
 LOD_PI_2 = params.p_LOD_PI_2; 
 
      %take out all Nans
    for i = 1:2
        index = isnan(data_PI_2(i,:));
        time_temp = time_PI_2(i,:);
        data_temp  = data_PI_2(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_PI_2(i).time = time_temp; 
        indiv_PI_2(i).data = data_temp; 
        indiv_PI_2(i).LOD =  LOD_PI_2(i,1);
        indiv_PI_2(i).c = zeros(1, length(indiv_PI_2(i).data)); 
             for j = 1:length(indiv_PI_2(i).data)
                        if data_temp(j) <= indiv_PI_2(i).LOD      
                             indiv_PI_2(i).c(j) = 0;  
                        else
                             indiv_PI_2(i).c(j)  = 1; 
                        end
            end
    end
 
 time_PI_3 = params.p_time_PI_3; 
 data_PI_3 = params.p_data_PI_3; 
 LOD_PI_3 = params.p_LOD_PI_3 ; 
 
  for i = 1:2
        index = isnan(data_PI_3(i,:));
        time_temp = time_PI_3(i,:);
        data_temp  = data_PI_3(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_PI_3(i).time = time_temp; 
        indiv_PI_3(i).data = data_temp; 
        indiv_PI_3(i).LOD =  LOD_PI_3(i,1);
        indiv_PI_3(i).c = zeros(1, length(indiv_PI_3(i).data)); 
             for j = 1:length(indiv_PI_3(i).data)
                        if data_temp(j) <= indiv_PI_3(i).LOD      
                             indiv_PI_3(i).c(j) = 0;  
                        else
                             indiv_PI_3(i).c(j)  = 1; 
                        end
            end
  end
 
 time_SI_DF_1 = params.p_time_SI_DF_1; 
 data_SI_DF_1 = params.p_data_SI_DF_1; 
 LOD_SI_DF_1 = params.p_LOD_SI_DF_1;
 
  for i = 1:29
        index = isnan(data_SI_DF_1(i,:));
        time_temp = time_SI_DF_1(i,:);
        data_temp  = data_SI_DF_1(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_SI_DF_1(i).time = time_temp; 
        indiv_SI_DF_1(i).data = data_temp; 
        indiv_SI_DF_1(i).LOD =  LOD_SI_DF_1(i,1);
        indiv_SI_DF_1(i).c = zeros(1, length(indiv_SI_DF_1(i).data)); 
             for j = 1:length(indiv_SI_DF_1(i).data)
                        if data_temp(j) <= indiv_SI_DF_1(i).LOD      
                             indiv_SI_DF_1(i).c(j) = 0;  
                        else
                             indiv_SI_DF_1(i).c(j)  = 1; 
                        end
            end
   end
 
 time_SI_DF_2 = params.p_time_SI_DF_2; 
 data_SI_DF_2 = params.p_data_SI_DF_2; 
 LOD_SI_DF_2 = params.p_LOD_SI_DF_2;
 
   for i = 1:3
        index = isnan(data_SI_DF_2(i,:));
        time_temp = time_SI_DF_2(i,:);
        data_temp  = data_SI_DF_2(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_SI_DF_2(i).time = time_temp; 
        indiv_SI_DF_2(i).data = data_temp; 
        indiv_SI_DF_2(i).LOD =  LOD_SI_DF_2(i,1);
        indiv_SI_DF_2(i).c = zeros(1, length(indiv_SI_DF_2(i).data)); 
             for j = 1:length(indiv_SI_DF_2(i).data)
                        if data_temp(j) <= indiv_SI_DF_2(i).LOD      
                             indiv_SI_DF_2(i).c(j) = 0;  
                        else
                             indiv_SI_DF_2(i).c(j)  = 1; 
                        end
            end
   end
 
 time_SI_DF_3 = params.p_time_SI_DF_3; 
 data_SI_DF_3 = params.p_data_SI_DF_3; 
 LOD_SI_DF_3 = params.p_LOD_SI_DF_3;
 
  for i = 1:2
        index = isnan(data_SI_DF_3(i,:));
        time_temp = time_SI_DF_3(i,:);
        data_temp  = data_SI_DF_3(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_SI_DF_3(i).time = time_temp; 
        indiv_SI_DF_3(i).data = data_temp; 
        indiv_SI_DF_3(i).LOD =  LOD_SI_DF_3(i,1);
        indiv_SI_DF_3(i).c = zeros(1, length(indiv_SI_DF_3(i).data)); 
             for j = 1:length(indiv_SI_DF_3(i).data)
                        if data_temp(j) <= indiv_SI_DF_3(i).LOD      
                             indiv_SI_DF_3(i).c(j) = 0;  
                        else
                             indiv_SI_DF_3(i).c(j)  = 1; 
                        end
            end
   end
 
 time_SI_DHF_1 = params.p_time_SI_DHF_1; 
 data_SI_DHF_1 = params.p_data_SI_DHF_1; 
 LOD_SI_DHF_1 = params.p_LOD_SI_DHF_1;
 
   for i = 1:10
        index = isnan(data_SI_DHF_1(i,:));
        time_temp = time_SI_DHF_1(i,:);
        data_temp  = data_SI_DHF_1(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_SI_DHF_1(i).time = time_temp; 
        indiv_SI_DHF_1(i).data = data_temp; 
        indiv_SI_DHF_1(i).LOD =  LOD_SI_DHF_1(i,1);
        indiv_SI_DHF_1(i).c = zeros(1, length(indiv_SI_DHF_1(i).data)); 
             for j = 1:length(indiv_SI_DHF_1(i).data)
                        if data_temp(j) <= indiv_SI_DHF_1(i).LOD      
                             indiv_SI_DHF_1(i).c(j) = 0;  
                        else
                             indiv_SI_DHF_1(i).c(j)  = 1; 
                        end
            end
   end
 
 time_SI_DHF_2 = params.p_time_SI_DHF_2; 
 data_SI_DHF_2 = params.p_data_SI_DHF_2; 
 LOD_SI_DHF_2 = params.p_LOD_SI_DHF_2;
 
    for i = 1:6
        index = isnan(data_SI_DHF_2(i,:));
        time_temp = time_SI_DHF_2(i,:);
        data_temp  = data_SI_DHF_2(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_SI_DHF_2(i).time = time_temp; 
        indiv_SI_DHF_2(i).data = data_temp; 
        indiv_SI_DHF_2(i).LOD =  LOD_SI_DHF_2(i,1);
        indiv_SI_DHF_2(i).c = zeros(1, length(indiv_SI_DHF_2(i).data)); 
             for j = 1:length(indiv_SI_DHF_2(i).data)
                        if data_temp(j) <= indiv_SI_DHF_2(i).LOD      
                             indiv_SI_DHF_2(i).c(j) = 0;  
                        else
                             indiv_SI_DHF_2(i).c(j)  = 1; 
                        end
            end
   end
 
 time_SI_DHF_3 = params.p_time_SI_DHF_3; 
 data_SI_DHF_3 = params.p_data_SI_DHF_3; 
 LOD_SI_DHF_3 = params.p_LOD_SI_DHF_3;
 
    for i = 1:2
        index = isnan(data_SI_DHF_3(i,:));
        time_temp = time_SI_DHF_3(i,:);
        data_temp  = data_SI_DHF_3(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_SI_DHF_3(i).time = time_temp; 
        indiv_SI_DHF_3(i).data = data_temp; 
        indiv_SI_DHF_3(i).LOD =  LOD_SI_DHF_3(i,1);
        indiv_SI_DHF_3(i).c = zeros(1, length(indiv_SI_DHF_3(i).data)); 
             for j = 1:length(indiv_SI_DHF_3(i).data)
                        if data_temp(j) <= indiv_SI_DHF_3(i).LOD      
                             indiv_SI_DHF_3(i).c(j) = 0;  
                        else
                             indiv_SI_DHF_3(i).c(j)  = 1; 
                        end
            end
   end

 accep_temp = 0; 
 a_prop = 0; 
 
 parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
 params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20, 1e-3, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0.15, 1];
 
 my_posterior = zeros(iterations, 1); 
 
 my_cov = diag(1e-3.*ones(size(temp_chain, 2),1));

   first_posterior = multi_level_posterior_peak(parameters, temp_chain(1,:), k, indiv_PI_1,...
    indiv_PI_2, indiv_PI_3, indiv_SI_DF_1, indiv_SI_DF_2, indiv_SI_DF_3, indiv_SI_DHF_1, indiv_SI_DHF_2, indiv_SI_DHF_3);

    my_posterior(1, 1) = first_posterior;

l = 50000; 
m = 150000; 
z = 30000; 

 for i = 2:iterations
     
     if i > l && i <= m 
         my_cov = 5e-2*cov(temp_chain(z:(i-1), :));
     end
     
     if i > m
          my_cov = 5e-2*cov(temp_chain(z:m, :));
     end
    
     
     if i <= l
        means = temp_chain(i-1, :)./sd;
        cur_proposal = proposal(means, my_cov).*sd;
 
        if k == 1 || k == 4 || k == 5
            temp_proposal = cur_proposal; 
            temp_proposal(7) = cur_proposal(1) + cur_proposal(7); 
            temp_proposal(8) = cur_proposal(1) + cur_proposal(8); 
            
             while min(temp_proposal) <= 0 
                 cur_proposal = proposal(means, my_cov).*sd;
                 temp_proposal = cur_proposal; 
                 temp_proposal(7) = cur_proposal(1) + cur_proposal(7); 
                 temp_proposal(8) = cur_proposal(1) + cur_proposal(8);
             end

        elseif k == 2
            temp_proposal = cur_proposal; 
            temp_proposal(7) = cur_proposal(3) + cur_proposal(7); 
            temp_proposal(8) = cur_proposal(3) + cur_proposal(8); 
            
             while min(temp_proposal) <= 0 
                 cur_proposal = proposal(means, my_cov).*sd;
                 temp_proposal = cur_proposal; 
                 temp_proposal(7) = cur_proposal(3) + cur_proposal(7); 
                 temp_proposal(8) = cur_proposal(3) + cur_proposal(8);
             end
             
        elseif k == 3
            temp_proposal = cur_proposal; 
            temp_proposal(7) = cur_proposal(6) + cur_proposal(7); 
            temp_proposal(8) = cur_proposal(6) + cur_proposal(8); 
            
             while min(temp_proposal) <= 0 
                 cur_proposal = proposal(means, my_cov).*sd;
                 temp_proposal = cur_proposal; 
                 temp_proposal(7) = cur_proposal(6) + cur_proposal(7); 
                 temp_proposal(8) = cur_proposal(6) + cur_proposal(8);
             end
        end
     else
         
         means = temp_chain(i-1, :);
        cur_proposal = proposal(means, my_cov);
        
        if k == 1 || k == 4 || k == 5
            temp_proposal = cur_proposal; 
            temp_proposal(7) = cur_proposal(1) + cur_proposal(7); 
            temp_proposal(8) = cur_proposal(1) + cur_proposal(8); 
            
             while min(temp_proposal) <= 0 
                 cur_proposal = proposal(means, my_cov);
                 temp_proposal = cur_proposal; 
                 temp_proposal(7) = cur_proposal(1) + cur_proposal(7); 
                 temp_proposal(8) = cur_proposal(1) + cur_proposal(8);
             end
        elseif k == 2
            temp_proposal = cur_proposal; 
            temp_proposal(7) = cur_proposal(3) + cur_proposal(7); 
            temp_proposal(8) = cur_proposal(3) + cur_proposal(8); 
            
             while min(temp_proposal) <= 0 
                 cur_proposal = proposal(means, my_cov);
                 temp_proposal = cur_proposal; 
                 temp_proposal(7) = cur_proposal(3) + cur_proposal(7); 
                 temp_proposal(8) = cur_proposal(3) + cur_proposal(8);
             end
        elseif k == 3
            temp_proposal = cur_proposal; 
            temp_proposal(7) = cur_proposal(6) + cur_proposal(7); 
            temp_proposal(8) = cur_proposal(6) + cur_proposal(8); 
            
             while min(temp_proposal) <= 0 
                 cur_proposal = proposal(means, my_cov);
                 temp_proposal = cur_proposal; 
                 temp_proposal(7) = cur_proposal(6) + cur_proposal(7); 
                 temp_proposal(8) = cur_proposal(6) + cur_proposal(8);
             end
        end
     end

    %the posterior 
    
    current = multi_level_posterior_peak(parameters, cur_proposal, k, indiv_PI_1, indiv_PI_2, indiv_PI_3, indiv_SI_DF_1, indiv_SI_DF_2, indiv_SI_DF_3, indiv_SI_DHF_1, indiv_SI_DHF_2, indiv_SI_DHF_3);
    temporary = my_posterior(i - 1, 1);
    
if k == 1 || k == 4 || k == 5
    temp_past_chain = temp_chain(i-1,:);
    temp_past_chain(7) = temp_past_chain(1) + temp_past_chain(7); 
    temp_past_chain(8) = temp_past_chain(1) + temp_past_chain(8); 
    cdf_current = log(mvncdf(temp_proposal));
    cdf_temporary = log(mvncdf(temp_past_chain));
elseif k == 2
    temp_past_chain = temp_chain(i-1,:);
    temp_past_chain(7) = temp_past_chain(3) + temp_past_chain(7); 
    temp_past_chain(8) = temp_past_chain(3) + temp_past_chain(8); 
    cdf_current = log(mvncdf(temp_proposal));
    cdf_temporary = log(mvncdf(temp_past_chain));
elseif k == 3
    temp_past_chain = temp_chain(i-1,:);
    temp_past_chain(7) = temp_past_chain(6) + temp_past_chain(7); 
    temp_past_chain(8) = temp_past_chain(6) + temp_past_chain(8); 
    cdf_current = log(mvncdf(temp_proposal));
    cdf_temporary = log(mvncdf(temp_past_chain));
end
         
    probab = exp(current + cdf_temporary - temporary - cdf_current); %including cdf because of truncated mean
    
     if rand(1) < probab
         temp_chain(i, :) = cur_proposal;
         my_posterior(i) = current;
          i
         if i > m
             a_prop = a_prop + 1;
         end
         accep_temp = accep_temp + 1; 
     else
         temp_chain(i, :) = temp_chain(i-1,:); 
         my_posterior(i) = temporary; 
     end 
     
     if mod(i, 5000) == 0
        if k == 1
            save('temp_multi_level_chain_1_peak', 'temp_chain')
            save('temp_multi_level_posterior_1_peak', 'my_posterior')
        elseif k == 2
            save('temp_multi_level_chain_2_peak', 'temp_chain')
            save('temp_multi_level_posterior_2_peak', 'my_posterior')
        elseif k == 3
            save('temp_multi_level_chain_3_peak', 'temp_chain')
            save('temp_multi_level_posterior_3_peak', 'my_posterior')
        elseif k == 4
            save('temp_multi_level_chain_4_peak', 'temp_chain')
            save('temp_multi_level_posterior_4_peak', 'my_posterior')
        elseif k == 5
            save('temp_multi_level_chain_5_peak', 'temp_chain')
            save('temp_multi_level_posterior_5_peak', 'my_posterior')
        end
     end
   
 end
  
 if iterations > m
    accep_prop = a_prop/(iterations - m)
 else
    accep_prop = accep_temp/iterations
 end
        
chain = temp_chain;

 end
 
 