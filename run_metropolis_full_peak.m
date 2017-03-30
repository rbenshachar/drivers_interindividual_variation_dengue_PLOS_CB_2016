 function [chain, accep_prop, my_posterior] = run_metropolis_full_peak(startvalue, iterations, params, k, sd)
 
     if k == 0
         temp_chain = zeros(iterations, 5); 
         temp_chain(1,:) = startvalue;
     elseif k == 1
         temp_chain = zeros(iterations, 6); 
         temp_chain(1,:) = startvalue;
     end

     time_PI = params.p_time_PI; 
     data_PI = params.p_data_PI; 
     LOD_PI = params.p_LOD_PI; 

     time_SI = params.p_time_SI; 
     data_SI = params.p_data_SI; 
     LOD_SI = params.p_LOD_SI;
     
     %take out all Nans
    for i = 1:size(data_PI, 1)
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

    for i = 1:size(data_SI, 1)
        index = isnan(data_SI(i,:));
        time_temp = time_SI(i,:);
        data_temp  = data_SI(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_SI(i).time = time_temp; 
        indiv_SI(i).data = data_temp; 
        indiv_SI(i).LOD =  LOD_SI(i,1);
        indiv_SI(i).c = zeros(1, length(indiv_SI(i).data)); 
            for j = 1:length(indiv_SI(i).data)
                        if data_temp(j) <= indiv_SI(i).LOD    
                             indiv_SI(i).c(j) = 0;  
                        else
                             indiv_SI(i).c(j)  = 1; 
                        end
            end
    end

     
         accep_temp = 0; 

         parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
         params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20, 1e-3, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0, 1];
     
         my_posterior = zeros(iterations, 1); 
         a_prop = 0; 
         my_cov = diag(5e-3.*ones(size(temp_chain, 2),1));

         %checking difference between model simulation T and data T
         
         
        first_posterior = posterior_peak(parameters,temp_chain(1,:), k, indiv_PI, indiv_SI);

        my_posterior(1, 1) = first_posterior;       
      
    l = 50000; 
    m = 150000; 
    z = 30000; 
        
 for i = 2:iterations 
     
     if i >= l && i <= m 
         my_cov = cov(temp_chain(z:(i-1), :));
     end
     
     if i > m
          my_cov = cov(temp_chain(z:m, :));
     end
    
     
     if i < l
     
        means = temp_chain(i-1, :)./sd;
        cur_proposal = proposal(means, my_cov).*sd;
   
           while min(cur_proposal) <= 0 
                 cur_proposal = proposal(means, my_cov).*sd;
           end
       
     else
         
         means = temp_chain(i-1, :);
        cur_proposal = proposal(means, my_cov);

         while min(cur_proposal) <= 0 
             cur_proposal = proposal(means, my_cov);
         end
     end
     

    %the posterior 
    [current_posterior] = posterior_peak(parameters, cur_proposal, k, indiv_PI, indiv_SI); 
    current =  current_posterior;
    
    temporary = my_posterior(i - 1, 1);

    cdf_current = log(mvncdf(cur_proposal));
    cdf_temporary = log(mvncdf(temp_chain(i-1,:))); 
         
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
         if k == 0
            save('temp_chain_0_newest', 'temp_chain')
            save('temp_posterior_0_newest', 'my_posterior')
            i
         elseif k == 1
            save('temp_chain_1_peak', 'temp_chain')
            save('temp_posterior_1_peak', 'my_posterior')
            i
         else
            save('temp_chain_2_new', 'temp_chain')
            save('temp_posterior_2_new', 'my_posterior')
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
 
 