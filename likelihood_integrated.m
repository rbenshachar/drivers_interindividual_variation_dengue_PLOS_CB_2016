function [log_L] = likelihood_integrated(parameters, x, step_size, indiv, n)

%#codegen

parameters(10) = x(1); %beta
parameters(11) = x(2); %kappa
parameters(12) = x(3); %q
IP_g = x(4); %IPg
parameters(13) = x(5); %deltaT
parameters(16) = x(6); %qT
parameters(18) = x(7); %sigma
parameters(5) = x(8); %Vinit 
parameters(19) = x(9); %sigma_e
 
lik_ind =  zeros(1,n);

[T, Y] = eulers_method(0.001, parameters);

check = min(Y(3,:));
   
V = zeros(1, size(Y, 2));  
Ttemp = V; 

if check > 0 
     V = log10(Y(3, :));
     a = 0.001;
    Ttemp = round((1/a).*T);
end

p =[0.001, 0.999];
crit = logninv(p,log(IP_g),parameters(18));

%simpson's rule
increment = (crit(2) - crit(1))/step_size; 
IP_range = crit(1):(increment/2):crit(2); 
f = lognpdf(IP_range, log(IP_g), parameters(18)); 

for j = 1:n %for each individual           
    data_temp = indiv(j).data;
    time_temp = indiv(j).time; 
    LOD =  indiv(j).LOD;

    y1 = new_likelihood_simpson(parameters, time_temp, data_temp, LOD, V, Ttemp, IP_range, f, indiv(j).c);
    lik_ind(j) = log(y1);
end

log_L = sum(lik_ind);

end