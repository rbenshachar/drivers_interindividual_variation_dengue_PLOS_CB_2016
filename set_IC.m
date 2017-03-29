function IC = set_IC(n, parameters)

%setting IC using Latin Hypercube Sampling for checking global stability of parameters

%n is number of initial conditions
%parameters is vector of parameters

num_p = length(parameters);
IC = zeros(num_p, n); 

%varying params plus/minus 400%
a = 0.5;  %lower bound
b = 1.5;  %upper bound   

%Non-replacement sampling
increment = a:((b-a)/(n-1)):b;

%x is increment for each parameter
x = zeros(num_p, n); 

for i = 1:num_p
    x(i,:) = randperm(n); 
    IC(i,:) = parameters(i).*increment(x(i,:)); 
end

