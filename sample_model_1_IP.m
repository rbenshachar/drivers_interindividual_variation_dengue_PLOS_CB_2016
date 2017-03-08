function void = sample_model_1_IP()

%sample individual IPs from model 1 posterior using inverse transform sampling
%save sampled IP distribution stratified by CM and immune status
%calculate the posterior and cdf for each individual IP sampled
%Used for Fig.2a,b

close all;
load('params')

load('chain_1_IC_1_final')
load('my_posterior_1_IC_1_final') 

l = 150000; n = 300000; 
my_chain = chain; 
params.IP = 5.9; 
params.deltaT = 1e-6; 

%use parameter values that resulted in maximum posterior value
temp_val  = find(my_posterior == max(my_posterior(l:n))); 
params.beta  = my_chain(temp_val(1), 1); 
params.kappa = my_chain(temp_val(1), 2); 
params.q  = my_chain(temp_val(1), 3);
params.sigma = my_chain(temp_val(1), 4);
params.Vinit = 10^(my_chain(temp_val(1), 6));

parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
 params.Ninit, params.alpha, params.omega, params.d, 0,0,0, params.deltaT,params.dT,params.Tinit];

%Setting endpoints for distribution from which to draw IP from
p =[ 0.001, 0.999];
crit = logninv(p,log(params.IP), params.sigma);
LB = round(crit(1), 2); 
UB = round(crit(2), 2); 
IP = LB:.01:UB; 

PI_posterior = zeros(26, length(IP)); 
data = params.data_PI; 
time = params.time_PI; 
temp_LOD = params.LOD_PI; 
cdf_PI = PI_posterior;
IP_val_PI = zeros(26, 1); 

parameters(10) = params.beta; %beta
parameters(11) = params.kappa; %kappa
parameters(12) = params.q; %q
parameters(13) = 0; %deltaT
parameters(16) = 0; %qT
parameters(5) = params.Vinit;

[T, Y] = eulers_method(0.001, parameters);

check = min(Y(3,:));
   
V = zeros(1, size(Y, 2));  
Ttemp = V; 

if check > 0 
     V = log10(Y(3, :));
     a = 0.001;
    Ttemp = round((1/a).*T);
end

for j = 1:size(data, 1) %for each individual           
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
    counter = 0;

    for i = 1:length(IP)
        PI_posterior(j,i) =  log(lognpdf(IP(i),log(5.9), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i));

        counter = counter +  exp(PI_posterior(j,i));
        cdf_PI(j,i) = counter;
    end
    
    temp = sum(exp(PI_posterior(j,:))); 
    cdf_PI(j,:) = cdf_PI(j,:)./temp; 
end

%sampling IP for each individual based on posterior given all possible IPs
for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_PI(j,:) <= val); 
    IP_val_PI(j) = IP(vec(end));
end

%SI DF and DHF
    parameters(13) = params.deltaT;
    parameters(16) = my_chain(temp_val(1), 5); %qT
    params.qT = parameters(16);

    [T, Y] = eulers_method(0.001, parameters);

    check = min(Y(3,:));

    V = zeros(1, size(Y, 2));  
    Ttemp = V; 

    if check > 0 
         V = log10(Y(3, :));
         a = 0.001;
        Ttemp = round((1/a).*T);
    end

data = params.data_SI_DF; 
time = params.time_SI_DF; 
temp_LOD = params.LOD_SI_DF; 

SI_DF_posterior = zeros(138, length(IP)); 
cdf_SI_DF = SI_DF_posterior;
IP_val_SI_DF= zeros(138, 1); 

for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
   counter = 0;

    for i = 1:length(IP)
        SI_DF_posterior(j,i) =  log(lognpdf(IP(i),log(5.9), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i), params);

        counter = counter +  exp(SI_DF_posterior(j,i));
        cdf_SI_DF(j,i) = counter;
    end
    
    temp = sum(exp(SI_DF_posterior(j,:))); 
    cdf_SI_DF(j,:) = cdf_SI_DF(j,:)./temp; 

end

%sampling
for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_SI_DF(j,:) <= val); 
    IP_val_SI_DF(j) = IP(vec(end));
end

data = params.data_SI_DHF; 
time = params.time_SI_DHF; 
temp_LOD = params.LOD_SI_DHF;

SI_DHF_posterior = zeros(64, length(IP)); 
cdf_SI_DHF = SI_DHF_posterior;
IP_val_SI_DHF= zeros(64, 1); 

for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
   counter = 0;

    for i = 1:length(IP)
        SI_DHF_posterior(j,i) =  log(lognpdf(IP(i),log(5.9), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i), params);

        counter = counter +  exp(SI_DHF_posterior(j,i));
        cdf_SI_DHF(j,i) = counter;
    end
    temp = sum(exp(SI_DHF_posterior(j,:))); 
    cdf_SI_DHF(j,:) = cdf_SI_DHF(j,:)./temp;
end

%sampling
for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_SI_DHF(j,:) <= val); 
    IP_val_SI_DHF(j) = IP(vec(end));
end

save('IP_val_PI_1', 'IP_val_PI');
save('IP_val_SI_DF_1', 'IP_val_SI_DF');
save('IP_val_SI_DHF_1', 'IP_val_SI_DHF');

