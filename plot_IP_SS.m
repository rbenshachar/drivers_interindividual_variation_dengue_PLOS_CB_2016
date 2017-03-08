function void = plot_IP_SS(k)

close all;

if k == 1 
    
    load('chain_SS_1_IC_1_final')
    load('my_posterior_SS_1_IC_1_final')
    l = 150000; n = 300000;  
elseif k == 2
    load('chain_SS_2_IC_1_final') 
    load('my_posterior_SS_2_IC_1_final')
    l = 150000; n = 300000;  
else
    load('chain_SS_3_IC_1_final') 
    load('my_posterior_SS_3_IC_1_final')
    l = 150000; n = 300000;  
end

my_chain = chain_multi_level;
my_posterior = my_posterior_multi_level;

load('params');

parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20, 1e-3, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0, 1];

temp  = find(my_posterior == max(my_posterior(l:n))); 

%param values
if k== 1 
    params.beta_1 = my_chain(temp(1), 1); 
    params.kappa  = my_chain(temp(1), 2); 
    params.q = my_chain(temp(1), 3); 
    params.Vinit = 10^(my_chain(temp(1), 8));
    params.sigma = my_chain(temp(1), 4);
    params.qT  = my_chain(temp(1), 5); 
    params.beta_2 = my_chain(temp(1), 1) + my_chain(temp(1), 6); 
    params.beta_3 = my_chain(temp(1), 1) + my_chain(temp(1), 7);
elseif k == 2
    params.beta = my_chain(temp(1), 1); 
    params.kappa  = my_chain(temp(1), 2); 
    params.q_1 = my_chain(temp(1), 3); 
    params.Vinit = 10^(my_chain(temp(1), 8));
    params.sigma = my_chain(temp(1), 4);
    params.qT  = my_chain(temp(1), 5); 
    params.q_2 = my_chain(temp(1), 3) + my_chain(temp(1), 6); 
    params.q_3 = my_chain(temp(1), 3) + my_chain(temp(1), 7);
else
    params.beta = my_chain(temp(1), 1); 
    params.kappa  = my_chain(temp(1), 2); 
    params.q = my_chain(temp(1), 3); 
    params.Vinit = 10^(my_chain(temp(1), 8));
    params.sigma = my_chain(temp(1), 4);
    params.qT_1  = my_chain(temp(1), 5); 
    params.qT_2 = my_chain(temp(1), 5) + my_chain(temp(1), 6); 
    params.qT_3 = my_chain(temp(1), 5) + my_chain(temp(1), 7);
end

%finding IP_js
p =[ 0.001, 0.999];
params.IP = 5.9; 
crit = logninv(p,log(params.IP), params.sigma);
LB = round(crit(1), 2); 
UB = round(crit(2), 2); 

IP = LB:.01:UB; 

%finding the IP_js

%DENV-1 PI
if k == 1
    params.beta = params.beta_1;
elseif k ==2
    params.q = params.q_1;
else
    params.qT = params.qT_1;
end

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

data = params.data_PI_1; 
time = params.time_PI_1; 
temp_LOD = params.LOD_PI_1; 

IP_val_PI_1 = zeros(15, 1); 
PI_posterior_1 = zeros(15, length(IP)); 
cdf_PI_1 = PI_posterior_1;

for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
    counter = 0;

    for i = 1:length(IP)
        PI_posterior_1(j,i) =  log(lognpdf(IP(i),log(5.9), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i));

        counter = counter +  exp(PI_posterior_1(j,i));
        cdf_PI_1(j,i) = counter;
    end
    
    temp = sum(exp(PI_posterior_1(j,:))); 
    cdf_PI_1(j,:) = cdf_PI_1(j,:)./temp;
end

%sampling
for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_PI_1(j,:) <= val); 
    IP_val_PI_1(j) = IP(vec(end)); 
end

%DENV-2 PI
if k == 1
    params.beta = params.beta_2;
elseif k ==2
    params.q = params.q_2;
else
    params.qT = params.qT_2;
end

parameters(10) = params.beta; %beta
parameters(12) = params.q; %q
parameters(13) = 0; %deltaT
parameters(16) = 0; %qT

[T, Y] = eulers_method(0.001, parameters);

check = min(Y(3,:));
   
V = zeros(1, size(Y, 2));  
Ttemp = V; 

if check > 0 
     V = log10(Y(3, :));
     a = 0.001;
    Ttemp = round((1/a).*T);
end

data = params.data_PI_2; 
time = params.time_PI_2; 
temp_LOD = params.LOD_PI_2; 

IP_val_PI_2 = zeros(5, 1); 
PI_posterior_2 = zeros(5, length(IP)); 
cdf_PI_2 = PI_posterior_2;

for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
   counter = 0;

    for i = 1:length(IP)
        PI_posterior_2(j,i) =  log(lognpdf(IP(i),log(5.9), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i));

        counter = counter +  exp(PI_posterior_2(j,i));
        cdf_PI_2(j,i) = counter;
    end
    
    temp = sum(exp(PI_posterior_2(j,:))); 
    cdf_PI_2(j,:) = cdf_PI_2(j,:)./temp;
end

%sampling
for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_PI_2(j,:) <= val); 
    IP_val_PI_2(j) = IP(vec(end));
end

%DENV-3 PI
if k == 1
    params.beta = params.beta_3;
elseif k ==2
    params.q = params.q_3;
else
    params.qT = params.qT_3;
end

parameters(10) = params.beta; %beta
parameters(11) = params.kappa; %kappa
parameters(12) = params.q; %q
parameters(13) = 0; %deltaT
parameters(16) = 0; %qT

[T, Y] = eulers_method(0.001, parameters);

check = min(Y(3,:));
   
V = zeros(1, size(Y, 2));  
Ttemp = V; 

if check > 0 
     V = log10(Y(3, :));
     a = 0.001;
    Ttemp = round((1/a).*T);
end

data = params.data_PI_3; 
time = params.time_PI_3; 
temp_LOD = params.LOD_PI_3; 

IP_val_PI_3 = zeros(6, 1); 
PI_posterior_3 = zeros(6, length(IP)); 
cdf_PI_3 = PI_posterior_3;

for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
        counter = 0;

    for i = 1:length(IP)
        PI_posterior_3(j,i) =  log(lognpdf(IP(i),log(5.9), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i));

        counter = counter +  exp(PI_posterior_3(j,i));
        cdf_PI_3(j,i) = counter;
    end
    
    temp = sum(exp(PI_posterior_3(j,:))); 
    cdf_PI_3(j,:) = cdf_PI_3(j,:)./temp;
end

%sampling
for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_PI_3(j,:) <= val); 
    IP_val_PI_3(j) = IP(vec(end));
end

%DENV-1 SI 
if k == 1
    params.beta = params.beta_1;
elseif k ==2
    params.q = params.q_1;
else 
    params.qT = params.qT_1;
end

parameters(10) = params.beta; %beta
parameters(11) = params.kappa; %kappa
parameters(12) = params.q; %q
parameters(13) = params.deltaT; %deltaT
parameters(16) = params.qT; %qT

[T, Y] = eulers_method(0.001, parameters);

check = min(Y(3,:));
   
V = zeros(1, size(Y, 2));  
Ttemp = V; 

if check > 0 
     V = log10(Y(3, :));
     a = 0.001;
    Ttemp = round((1/a).*T);
end

data = params.data_SI_1; 
time = params.time_SI_1; 
temp_LOD = params.LOD_SI_1; 

IP_val_SI_1 = zeros(size(data, 1), 1); 
SI_posterior_1  = zeros(size(data, 1), length(IP)); 
cdf_SI_1 = SI_posterior_1;

for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
        counter = 0;

    for i = 1:length(IP)
        SI_posterior_1(j,i) =  log(lognpdf(IP(i),log(5.9), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i));

        counter = counter +  exp(SI_posterior_1(j,i));
        cdf_SI_1(j,i) = counter;
    end
    
    temp = sum(exp(SI_posterior_1(j,:))); 
    cdf_SI_1(j,:) = cdf_SI_1(j,:)./temp;
end

%sampling
for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_SI_1(j,:) <= val); 
    IP_val_SI_1(j) = IP(vec(end));

end

%DENV-2 SI 
if k == 1
    params.beta = params.beta_2;
elseif k == 2
    params.q = params.q_2;
else
    params.qT = params.qT_2;
end

parameters(10) = params.beta; %beta
parameters(11) = params.kappa; %kappa
parameters(12) = params.q; %q
parameters(13) = params.deltaT; %deltaT
parameters(16) = params.qT; %qT

[T, Y] = eulers_method(0.001, parameters);

check = min(Y(3,:));
   
V = zeros(1, size(Y, 2));  
Ttemp = V; 

if check > 0 
     V = log10(Y(3, :));
     a = 0.001;
    Ttemp = round((1/a).*T);
end

data = params.data_SI_2; 
time = params.time_SI_2; 
temp_LOD = params.LOD_SI_2; 

IP_val_SI_2 = zeros(size(data, 1), 1); 
SI_posterior_2 = zeros(length(IP), size(data, 1)); 
cdf_SI_2 = SI_posterior_2; 

for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
        counter = 0;

    for i = 1:length(IP)
        SI_posterior_2(j,i) =  log(lognpdf(IP(i),log(5.9), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i));

        counter = counter +  exp(SI_posterior_2(j,i));
        cdf_SI_2(j,i) = counter;
    end
    
    temp = sum(exp(SI_posterior_2(j,:))); 
    cdf_SI_2(j,:) = cdf_SI_2(j,:)./temp;
end

%sampling
for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_SI_2(j,:) <= val); 
    IP_val_SI_2(j) = IP(vec(end));
end

%DENV-3 SI
if k == 1
    params.beta = params.beta_3;
elseif k == 2
    params.q = params.q_3;
else
    params.qT = params.qT_3;
end

parameters(10) = params.beta; %beta
parameters(11) = params.kappa; %kappa
parameters(12) = params.q; %q
parameters(13) = params.deltaT; %deltaT
parameters(16) = params.qT; %qT

[T, Y] = eulers_method(0.001, parameters);

check = min(Y(3,:));
   
V = zeros(1, size(Y, 2));  
Ttemp = V; 

if check > 0 
     V = log10(Y(3, :));
     a = 0.001;
    Ttemp = round((1/a).*T);
end

data = params.data_SI_3; 
time = params.time_SI_3; 
temp_LOD = params.LOD_SI_3; 

IP_val_SI_3 = zeros(size(data,1), 1); 
SI_posterior_3 = zeros(length(IP), size(data,1)); 
cdf_SI_3 = SI_posterior_3; 

for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
         counter = 0;

    for i = 1:length(IP)
        SI_posterior_3(j,i) =  log(lognpdf(IP(i),log(5.9), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i));

        counter = counter +  exp(SI_posterior_3(j,i));
        cdf_SI_3(j,i) = counter;
    end
    
    temp = sum(exp(SI_posterior_3(j,:))); 
    cdf_SI_3(j,:) = cdf_SI_3(j,:)./temp;
end

%sampling
for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_SI_3(j,:) <= val); 
    IP_val_SI_3(j) = IP(vec(end));
end

if k == 1
    save('IP_val_PI_1_beta', 'IP_val_PI_1')
    save('IP_val_PI_2_beta', 'IP_val_PI_2')
    save('IP_val_PI_3_beta', 'IP_val_PI_3')
    save('IP_val_SI_1_beta', 'IP_val_SI_1')
    save('IP_val_SI_2_beta', 'IP_val_SI_2')
    save('IP_val_SI_3_beta', 'IP_val_SI_3')
elseif k == 2
    save('IP_val_PI_1_q', 'IP_val_PI_1')
    save('IP_val_PI_2_q', 'IP_val_PI_2')
    save('IP_val_PI_3_q', 'IP_val_PI_3')
    save('IP_val_SI_1_q', 'IP_val_SI_1')
    save('IP_val_SI_2_q', 'IP_val_SI_2')
    save('IP_val_SI_3_q', 'IP_val_SI_3')
elseif k == 3
    save('IP_val_PI_1_qT', 'IP_val_PI_1')
    save('IP_val_PI_2_qT', 'IP_val_PI_2')
    save('IP_val_PI_3_qT', 'IP_val_PI_3')
    save('IP_val_SI_1_qT', 'IP_val_SI_1')
    save('IP_val_SI_2_qT', 'IP_val_SI_2')
    save('IP_val_SI_3_qT', 'IP_val_SI_3')
end


%plot_IP_distribution
IP_val_PI = cat(1, IP_val_PI_1, IP_val_PI_2, IP_val_PI_3); 
IP_val_SI = cat(1, IP_val_SI_1, IP_val_SI_2, IP_val_SI_3); 

subplot(3,1,1); 
h = histogram(IP_val_PI); 
h.Normalization = 'probability';
h.NumBins = 20;
title('PI DF')

subplot(3,1,2)
h = histogram(IP_val_SI); 
h.Normalization = 'probability'; 
h.NumBins = 20;
title('SI DHF')
xlabel('individual IP')

subplot(3,1,3)
IP_val_PI_sort = sort(IP_val_PI);
plot(IP_val_PI_sort,(1:length(IP_val_PI_sort))./length(IP_val_PI_sort), 'm')
hold on; 
IP_val_SI_sort = sort(IP_val_SI);
plot(IP_val_SI_sort,(1:length(IP_val_SI_sort))./length(IP_val_SI_sort), 'b')
legend('PI', 'SI')
xlabel('individual IP')
ylabel('cumulative number of individuals with IP by group')

