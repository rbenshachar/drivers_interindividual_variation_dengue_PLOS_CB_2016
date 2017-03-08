function value = find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, x)

%returns the log-likelihood value for an individual with a specified IP

%INPUTS: 
%parameters
%time_temp: individual time-points of hospital admission
%data_temp: individual data
%LOD: LOD data
%V: model-simulated viremia
%Ttemp: simulated time vector
%x: individual IP value

%OUTPUT: Log-likelihood value

%vector if viral load measurements are above or below LOD - for easy calculation of likelihood
   c = zeros(1, length(data_temp)); 
   for i = 1:length(data_temp) %individual data points       
        if data_temp(i) <= LOD      
             c(i) = 0;  
        else
             c(i) = 1; 
        end
   end
  
 parameters(17) = x;
 parameters(19) = 1; 
  
 simulation_current = model_fit(parameters, time_temp, V, Ttemp);
  
 temp = get_likelihood(data_temp, LOD, simulation_current, parameters(19), c, 1-c); %likelihood not on log scale
 
 value = log(temp); 
