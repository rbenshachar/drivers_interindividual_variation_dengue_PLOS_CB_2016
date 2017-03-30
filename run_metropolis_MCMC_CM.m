 function [chain, accep_prop, my_posterior] = run_metropolis_MCMC_CM(startvalue, iterations, params, k, sd) %, current_val) %, current_val)
  
 if k == 1 %OAS
     temp_chain = zeros(iterations, 7); 
     temp_chain(1,:) = startvalue;
 elseif k  == 2 %ADE
     temp_chain = zeros(iterations, 7); 
     temp_chain(1,:) = startvalue;
 elseif k == 3 %alternative OAS model
     temp_chain = zeros(iterations, 8); 
     temp_chain(1,:) = startvalue;
 end
 
 num_params = size(temp_chain, 2);
  
 time_PI = params.time_PI; 
 data_PI = params.data_PI; 
 LOD_PI = params.LOD_PI; 
 
      %take out all Nans
    for i = 1:26
        index = isnan(data_PI(i,:));
        time_temp = time_PI(i,:);
        data_temp  = data_PI(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_PI(i).time = time_temp; 
        indiv_PI(i).data = data_temp; 
        indiv_PI(i).LOD =  LOD_PI(i,1);
        indiv_PI(i).c = zeros(1, length(indiv_PI(i).data)); 
             for j = 1:length(indiv_PI(i).data)
                        if data_temp(j) <= indiv_PI(i).LOD      
                             indiv_PI(i).c(j) = 0;  
                        else
                             indiv_PI(i).c(j)  = 1; 
                        end
            end
    end
    
   time_SI_DF = params.time_SI_DF; 
   data_SI_DF = params.data_SI_DF; 
   LOD_SI_DF = params.LOD_SI_DF;

    for i = 1:138
        index = isnan(data_SI_DF(i,:));
        time_temp = time_SI_DF(i,:);
        data_temp  = data_SI_DF(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_SI_DF(i).time = time_temp; 
        indiv_SI_DF(i).data = data_temp; 
        indiv_SI_DF(i).LOD = LOD_SI_DF(i,1);
        indiv_SI_DF(i).c = zeros(1, length(indiv_SI_DF(i).data)); 
            for j = 1:length(indiv_SI_DF(i).data)
                        if data_temp(j) <= indiv_SI_DF(i).LOD     
                             indiv_SI_DF(i).c(j) = 0;  
                        else
                             indiv_SI_DF(i).c(j)  = 1; 
                        end
            end
    end
 
 time_SI_DHF = params.time_SI_DHF; 
 data_SI_DHF = params.data_SI_DHF; 
 LOD_SI_DHF = params.LOD_SI_DHF;
    
   for i = 1:64
        index = isnan(data_SI_DHF(i,:));
        time_temp = time_SI_DHF(i,:);
        data_temp  = data_SI_DHF(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_SI_DHF(i).time = time_temp; 
        indiv_SI_DHF(i).data = data_temp; 
        indiv_SI_DHF(i).LOD =  LOD_SI_DHF(i,1);
        indiv_SI_DHF(i).c = zeros(1, length(indiv_SI_DHF(i).data)); 
           for j = 1:length(indiv_SI_DHF(i).data)
                        if data_temp(j) <= indiv_SI_DHF(i).LOD     
                             indiv_SI_DHF(i).c(j) = 0;  
                        else
                             indiv_SI_DHF(i).c(j)  = 1; 
                        end
            end
    end
 
 accep_temp = 0; 
 a_prop = 0; 
 
 my_posterior = zeros(iterations, 1); 

  parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
  params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20, 1e-3, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0.15, 1];

   my_cov = diag(1e-3.*ones(size(temp_chain, 2),1));

   first_posterior = posterior_CM(parameters, temp_chain(1,:), k, indiv_PI, indiv_SI_DF, indiv_SI_DHF);

   my_posterior(1, 1) = first_posterior; 

 l = 50000; 
 m = 150000; 
 z = 30000; 
 
 for i = 2:iterations
     
      if i >= l && i <= m 
         my_cov = cov(temp_chain(z:(i- 1), :));
     end
     
     if i > m
          my_cov = cov(temp_chain(z:m, :));
     end
    
     
     if i < l
     
        means = temp_chain(i-1, :)./sd;
        cur_proposal = proposal(means, my_cov).*sd;
   
           while min(cur_proposal(1:(num_params - 1))) <= 0 
                 cur_proposal = proposal(means, my_cov).*sd;
           end
       
     else
         
         means = temp_chain(i-1, :);
        cur_proposal = proposal(means, my_cov);
        

         while min(cur_proposal(1:(num_params - 1))) <= 0 
             cur_proposal = proposal(means, my_cov);
         end
     end

%the posterior 
    
current =  posterior_CM(parameters, cur_proposal, k, indiv_PI, indiv_SI_DF, indiv_SI_DHF);
temporary = my_posterior(i - 1, 1);

cdf_current = log(mvncdf(cur_proposal(1:(num_params - 1))));
cdf_temporary = log(mvncdf(temp_chain(i-1,(1:(num_params - 1)))));

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
            save('temp_chain_CM_1', 'temp_chain')
            save('temp_posterior_CM_1', 'my_posterior')
            i
         elseif k == 2
            save('temp_chain_CM_2', 'temp_chain')
            save('temp_posterior_CM_2', 'my_posterior')
            i
         elseif k == 3
            save('temp_chain_CM_3', 'temp_chain')
            save('temp_posterior_CM_3', 'my_posterior')
            i
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
 
 