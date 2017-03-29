 function [chain, accep_prop, my_posterior] = run_metropolis_full(startvalue, iterations, params, k, sd, p, s)

%startvalue is vector of parameter initial conditions
%iterations is number of iterations
%k is model identifier
%sd is standard deviations of initial parameter estimates
%p = number of primary infections
%s = number of secondary infection


     if k == 0 || k == 5
         temp_chain = zeros(iterations, 5);
         temp_chain(1,:) = startvalue;
     elseif k == 1 || k == 2 || k == 3
         temp_chain = zeros(iterations, 6); 
         temp_chain(1,:) = startvalue;
     elseif k == 4 
         temp_chain = zeros(iterations, 7); 
         temp_chain(1,:) = startvalue; %chapter 3 
     end

accep_temp = 0;

%initial posterior matrix
my_posterior = zeros(iterations, 1);
%initialize acceptance ratio
a_prop = 0;
%initial covariance matrix for proposal
my_cov = diag(5e-3.*ones(size(temp_chain, 2),1));

parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
         params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20, 1e-3, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0, 1];
    
 time_PI = params.time_PI; 
 data_PI = params.data_PI; 
 LOD_PI = params.LOD_PI;      
     
    for i = 1:p %26 
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
                             indiv_PI(i).c(j) = 1; 
                        end
            end
    end
     
 time_SI = params.time_SI; 
 data_SI = params.data_SI; 
 LOD_SI = params.LOD_SI;      
     
    for i = 1:s %202 
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
                             indiv_SI(i).c(j) = 1; 
                        end
            end
    end

 %initialize first posterior 
 first_posterior = posterior(parameters,temp_chain(1,:), k, indiv_PI, indiv_SI);

 my_posterior(1, 1) = first_posterior;
      
    l = 50000; 
    m = 150000; 
    z = 30000; 
        
 for i = 2:iterations 
     
     if i >= l && i <= m 
         if k == 4 || k == 5
            my_cov = 5e-2.*cov(temp_chain(z:(i-1), :));
         else
            my_cov = cov(temp_chain(z:(i-1), :)); %start adapting covariance matrix
         end
     end
     
     if i > m
         if k == 4 || k == 5
             my_cov = 5e-2.*cov(temp_chain(z:m, :));
         else
            my_cov = cov(temp_chain(z:m, :)); %fix covariance
         end
     end
    
     
     if i < l %start with small sd
     
        means = temp_chain(i-1, :)./sd;
        cur_proposal = proposal(means, my_cov).*sd;
   
        if k == 0
            while min(cur_proposal(1:4)) <= 0 
                 cur_proposal = proposal(means, my_cov).*sd;
           end
        else
           while min(cur_proposal(1:5)) <= 0 
                 cur_proposal = proposal(means, my_cov).*sd;
           end
        end
            
       
     else
         
         means = temp_chain(i-1, :);
        cur_proposal = proposal(means, my_cov);

        if k == 0
            while min(cur_proposal(1:4)) <= 0 
                 cur_proposal = proposal(means, my_cov);
             end
        else
             while min(cur_proposal(1:5)) <= 0 
                 cur_proposal = proposal(means, my_cov);
             end
        end
     end
     

    %the posterior 
    [current_posterior] = posterior(parameters, cur_proposal, k, indiv_PI, indiv_SI); 
    current =  current_posterior;
    
    temporary = my_posterior(i - 1, 1);

    %the cdfs are used to account for truncating proposal at 0
    if k == 0
        cdf_current = log(mvncdf(cur_proposal(1:4)));
        cdf_temporary = log(mvncdf(temp_chain(i-1,1:4))); 
    else
        cdf_current = log(mvncdf(cur_proposal(1:5)));
        cdf_temporary = log(mvncdf(temp_chain(i-1,1:5))); 
    end
         
    probab = exp(current + cdf_temporary - temporary - cdf_current); %including cdf because of truncated mean

    %accept or refect parameter sets
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
            save('temp_chain_0', 'temp_chain')
            save('temp_posterior_0', 'my_posterior')
            i
         elseif k == 1
            save('temp_chain_1', 'temp_chain')
            save('temp_posterior_1', 'my_posterior')
            i
         elseif k == 2
            save('temp_chain_2', 'temp_chain')
            save('temp_posterior_2', 'my_posterior')
            i
         elseif k == 3
              save('temp_chain_3', 'temp_chain')
            save('temp_posterior_3', 'my_posterior')
            i
         elseif k == 4
             save('temp_chain_4', 'temp_chain')
            save('temp_posterior_4', 'my_posterior')
            i 
       elseif k == 5
             save('temp_chain_5', 'temp_chain')
            save('temp_posterior_5', 'my_posterior')
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
 
 