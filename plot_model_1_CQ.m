%plot_model_1_CQ

%Figure S4 in PLOS CB paper
%viral load fits of model 1 fit to patients who did and did not receive
%Chloroquine separately

close all;
clear all; 

fig = figure;

load('params');

colors = colormap(cbrewer('qual', 'Set3', 5));

load('chain_1_IC_1_final_P')
load('my_posterior_1_IC_1_final_P');

l = 150000; m = 100; n = 300000; 
my_chain = chain; 

k = 100; 
samples = round(l + ((n-l))*rand(k, 1));

params.time_end = 15; 
     
viremia_PI_P = zeros(k, length(params.time_start:.01:params.time_end)); 
viremia_SI_P = zeros(k, length(params.time_start:.01:params.time_end));


 x = 1; y = 2;

for i = 1:length(samples)
    params.beta = my_chain(samples(i), 1); 
    params.kappa = my_chain(samples(i), 2); 
    params.q = my_chain(samples(i), 3); 
    params.sigma = my_chain(samples(i), 4); 
    params.qT = my_chain(samples(i), 5); 
    params.Vinit = 10^(my_chain(samples(i), 6)); 

    params.time_end = 15; 
    
     [T, Y] = ode45(@(t,y)PI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);  
     [T2, Y2] = ode45(@(t,y)SI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);     
    
     T_PI = T;
     viremia_PI_P(i,:) = log10(Y(:,3)); 
     
     T_SI = T2; 
     viremia_SI_P(i,:) = log10(Y2(:,3)); 
end

load('chain_1_IC_1_final_CQ')
load('my_posterior_1_IC_1_final_CQ');
my_chain = chain; 


viremia_PI_C = zeros(k, length(params.time_start:.01:params.time_end)); 
viremia_SI_C = zeros(k, length(params.time_start:.01:params.time_end));

for i = 1:length(samples)
    params.beta = my_chain(samples(i), 1); 
    params.kappa = my_chain(samples(i), 2); 
    params.q = my_chain(samples(i), 3); 
    params.sigma = my_chain(samples(i), 4); 
    params.qT = my_chain(samples(i), 5); 
    params.Vinit = 10^(my_chain(samples(i), 6)); 

    params.time_end = 15; 
    
     [T, Y] = ode45(@(t,y)PI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);  
     [T2, Y2] = ode45(@(t,y)SI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);     
    
     viremia_PI_C(i,:) = log10(Y(:,3)); 
     viremia_SI_C(i,:) = log10(Y2(:,3)); 
end

viremia_PI_P_low = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_PI_P_med = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_PI_P_high = zeros(length(params.time_start:.01:params.time_end), 1)';

viremia_SI_P_low = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_SI_P_med = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_SI_P_high = zeros(length(params.time_start:.01:params.time_end), 1)';

viremia_PI_CQ_low = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_PI_CQ_med = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_PI_CQ_high = zeros(length(params.time_start:.01:params.time_end), 1)';

viremia_SI_CQ_low = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_SI_CQ_med = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_SI_CQ_high = zeros(length(params.time_start:.01:params.time_end), 1)';


for j = 1:length(params.time_start:.01:params.time_end)
  
   viremia_PI_P_low(j) =  quantile(viremia_PI_P(:,j), 0.0250); 
   viremia_PI_P_med(j) = quantile(viremia_PI_P(:,j), 0.5); 
   viremia_PI_P_high(j) = quantile(viremia_PI_P(:,j), 0.975); 
   
   viremia_SI_P_low(j) =  quantile(viremia_SI_P(:,j), 0.0250); 
   viremia_SI_P_med(j) = quantile(viremia_SI_P(:,j), 0.5); 
   viremia_SI_P_high(j) = quantile(viremia_SI_P(:,j), 0.975); 
   
   viremia_PI_CQ_low(j) =  quantile(viremia_PI_C(:,j), 0.0250); 
   viremia_PI_CQ_med(j) = quantile(viremia_PI_C(:,j), 0.5); 
   viremia_PI_CQ_high(j) = quantile(viremia_PI_C(:,j), 0.975); 
   
   viremia_SI_CQ_low(j) =  quantile(viremia_SI_C(:,j), 0.0250); 
   viremia_SI_CQ_med(j) = quantile(viremia_SI_C(:,j), 0.5); 
   viremia_SI_CQ_high(j) = quantile(viremia_SI_C(:,j), 0.975); 
   
end

subplot(x,y,1)  
shadedplot(T_PI, viremia_PI_P_low,  viremia_PI_P_high, colors(4,:), 'none');
hold on;
shadedplot(T_PI, viremia_PI_CQ_low,  viremia_PI_CQ_high, colors(2,:), 'none');
hold on; 
plot(T_PI, viremia_PI_P_med, 'k', 'LineWidth', 2); 
hold on; 
plot(T_PI, viremia_PI_CQ_med, 'k-.', 'LineWidth', 2); 
text(0.1,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
xlabel('time since infection (days)')
ylabel('viral load (log copies/ml)')
%set(gca, 'FontSize', 16)
ylim([-4, 12])
    
subplot(x,y,2)  
shadedplot(T_SI, viremia_SI_P_low,  viremia_SI_P_high, colors(4,:), 'none');
hold on;
shadedplot(T_SI, viremia_SI_CQ_low,  viremia_SI_CQ_high, colors(2,:), 'none');
hold on; 
plot(T_SI, viremia_SI_P_med, 'k', 'LineWidth', 2); 
hold on; 
plot(T_SI, viremia_SI_CQ_med, 'k-.', 'LineWidth', 2); 
text(0.1,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
xlabel('time since infection (days)')
ylabel('viral load (log copies/ml)')
%set(gca, 'FontSize', 16)
ylim([-4, 12])

   fig.Units = 'inches';
   fig.PaperPosition = [1 1 7.5 3.5];
   fig.Position = fig.PaperPosition;

   print(fig,'figS4.tiff','-dtiff','-r300')


