% plot_best_SS_models

% FIGURE 6 and S6 IN PLOS CB PAPER

close all;
clear all; 

load('params');
params.time_end = 15; 

fig = figure;
l = 150000; m = 1; n = 300000;  
k = 100; 
samples = round(l + ((n-l))*rand(k, 1));

colors = colormap(cbrewer('qual', 'Set3', 5));

%model 1
load('chain_SS_1_IC_1_final')
load('my_posterior_SS_1_IC_1_final')

my_chain = chain_multi_level;
my_posterior = my_posterior_multi_level; 

temp  = find(my_posterior == max(my_posterior(l:n))); 

%  %IPs
%  %saved from plot_IP_SS(1)
load('IP_val_PI_1_beta')
load('IP_val_PI_2_beta')
load('IP_val_PI_3_beta')
load('IP_val_SI_1_beta')
load('IP_val_SI_2_beta')
load('IP_val_SI_3_beta')
 
%subplots
x = 2; y = 3; 

subplot(x,y,1)
for i = 1:size(params.time_PI_1, 1)
   plot(IP_val_PI_1(i) + params.time_PI_1(i,:),  params.data_PI_1(i,:),  'color', [0.5 0.5 0.5], 'LineWidth', 2);
    hold on; 
end
plot(0:15, 3.1761.*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','-.'); 
hold on; 
plot(0:15,4.1761.*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','-.');
text(0.1,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
%set(gca, 'FontSize', 16);
ylabel('viral load (log copies/ml)')
xlabel('time since infection (days)')
xlim([0, 15]); ylim([-5, 12]);

subplot(x,y,2)
for i = 1:size(params.time_PI_2, 1)
    plot(IP_val_PI_2(i) +params.time_PI_2(i,:),  params.data_PI_2(i,:),  'color', [0.5 0.5 0.5], 'LineWidth', 2);
    hold on; 
end
hold on; 
plot(0:15, 3.1761.*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','-.'); 
hold on; 
plot(0:15, 4.1761.*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','-.'); 
text(0.1,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
% set(gca, 'FontSize', 16);
ylabel('viral load (log copies/ml)')
xlabel('time since infection (days)')
xlim([0, 15]); ylim([-5, 12]);

subplot(x,y,3)
for i = 1:size(params.time_PI_3, 1)
    plot(IP_val_PI_3(i) +params.time_PI_3(i,:),  params.data_PI_3(i,:),  'color', [0.5 0.5 0.5], 'LineWidth', 2);
    hold on; 
end
plot(0:15,3.1761.*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','-.'); 
hold on; 
plot(0:15,4.1761.*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','-.');
text(0.1,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
%set(gca, 'FontSize', 16);
ylabel('viral load (log copies/ml)')
xlabel('time since infection (days)')
xlim([0, 15]); ylim([-5, 12]);

subplot(x,y,4)
for i = 1:size(params.time_SI_1, 1)
    plot(IP_val_SI_1(i) + params.time_SI_1(i,:), params.data_SI_1(i,:),  'color', [0.5 0.5 0.5], 'LineStyle','-.', 'LineWidth', 2);
    hold on; 
end
plot(0:15,3.1761.*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','-.'); 
hold on; 
plot(0:15,4.1761.*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','-.');
ylabel('viral load (log copies/ml)')
xlabel('time since infection (days)')
text(0.1,0.9,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
xlim([0, 15]); ylim([-5, 12]);
%set(gca, 'FontSize', 16); 


subplot(x,y,5)
for i = 1:size(params.time_SI_2, 1)
   plot(IP_val_SI_2(i) + params.time_SI_2(i,:),  params.data_SI_2(i,:),  'color', [0.5 0.5 0.5], 'LineStyle','-.', 'LineWidth', 2);
    hold on; 
end
plot(0:15,3.1761.*ones(16,1),'color', [0.5 0.5 0.5], 'LineStyle', '-.'); 
hold on; 
plot(0:15,4.1761.*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','-.');
xlabel('time since infection (days)')
ylabel('viral load (log copies/ml)')
text(0.1,0.9,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
xlim([0, 15]); ylim([-5, 12]);
%set(gca, 'FontSize', 16); 


subplot(x,y,6)
for i = 1:size(params.time_SI_3, 1)
    plot(IP_val_SI_3(i) + params.time_SI_3(i,:), params.data_SI_3(i,:),  'color', [0.5 0.5 0.5], 'LineStyle','-.', 'LineWidth', 2);
    hold on; 
end
plot(0:15,3.1761.*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','-.'); 
hold on; 
plot(0:15,4.1761.*ones(16,1),'color', [0.5 0.5 0.5], 'LineStyle', '-.');
xlabel('time since infection (days)')
ylabel('viral load (log copies/ml)')
text(0.1,0.9,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
xlim([0, 15]); ylim([-5, 12]);
%set(gca, 'FontSize', 16); 

viremia_1_PI = zeros(k, length(params.time_start:.1:params.time_end));
viremia_2_PI = zeros(k, length(params.time_start:.1:params.time_end));
viremia_3_PI = zeros(k, length(params.time_start:.1:params.time_end));

viremia_1_SI = zeros(k, length(params.time_start:.1:params.time_end));
viremia_2_SI = zeros(k, length(params.time_start:.1:params.time_end));
viremia_3_SI = zeros(k, length(params.time_start:.1:params.time_end));

T_1_SI = zeros(k, length(params.time_start:.1:params.time_end));
T_2_SI = zeros(k, length(params.time_start:.1:params.time_end));
T_3_SI = zeros(k, length(params.time_start:.1:params.time_end));

params.IP = 5.9; 

for i = 1:length(samples)
    params.beta_1 = my_chain(samples(i), 1); 
    params.kappa  = my_chain(samples(i), 2); 
    params.q = my_chain(samples(i), 3); 
    params.Vinit = 10^(my_chain(samples(i), 8));
    params.sigma = my_chain(samples(i), 4);
    params.qT  = my_chain(samples(i), 5); 
    params.beta_2 = my_chain(samples(i), 1) + my_chain(samples(i), 6); 
    params.beta_3 = my_chain(samples(i), 1) + my_chain(samples(i), 7);
    
   params.beta = params.beta_1;
   [TP1, YP1] = ode45(@(t,y)PI(t, y, params), params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);   
   [TS1, YS1] = ode45(@(t,y)SI(t, y, params), params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);        

   viremia_1_PI(i,:) = log10(YP1(:,3)); 
   viremia_1_SI(i,:) = log10(YS1(:,3)); 
   
   T_1_SI(i,:) = YS1(:,5); 
   
   params.beta = params.beta_2; 
   [TP2, YP2] = ode45(@(t,y)PI(t, y, params), params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);   
   [TS2, YS2] = ode45(@(t,y)SI(t, y, params), params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit ]);        

   viremia_2_PI(i,:) = log10(YP2(:,3)); 
   viremia_2_SI(i,:) = log10(YS2(:,3)); 
   
   T_2_SI(i,:) = YS2(:,5); 
   
   params.beta = params.beta_3; 
   [TP3, YP3] = ode45(@(t,y)PI(t, y, params), params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit ]);   
   [TS3, YS3] = ode45(@(t,y)SI(t, y, params), params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);        

   viremia_3_PI(i,:) = log10(YP3(:,3)); 
   viremia_3_SI(i,:) = log10(YS3(:,3)); 
   
   T_3_SI(i,:) = YS3(:,5); 
end

viremia_PI_1_low = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_PI_1_med = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_PI_1_high = zeros(length(params.time_start:.1:params.time_end), 1)';

viremia_SI_1_low = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_SI_1_med = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_SI_1_high = zeros(length(params.time_start:.1:params.time_end), 1)';

viremia_PI_2_low = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_PI_2_med = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_PI_2_high = zeros(length(params.time_start:.1:params.time_end), 1)';

viremia_SI_2_low = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_SI_2_med = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_SI_2_high = zeros(length(params.time_start:.1:params.time_end), 1)';

viremia_PI_3_low = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_PI_3_med = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_PI_3_high = zeros(length(params.time_start:.1:params.time_end), 1)';

viremia_SI_3_low = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_SI_3_med = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_SI_3_high = zeros(length(params.time_start:.1:params.time_end), 1)';

T_SI_1_low = zeros(length(params.time_start:.1:params.time_end), 1)';
T_SI_1_med = zeros(length(params.time_start:.1:params.time_end), 1)';
T_SI_1_high = zeros(length(params.time_start:.1:params.time_end), 1)';

T_SI_2_low = zeros(length(params.time_start:.1:params.time_end), 1)';
T_SI_2_med = zeros(length(params.time_start:.1:params.time_end), 1)';
T_SI_2_high = zeros(length(params.time_start:.1:params.time_end), 1)';

T_SI_3_low = zeros(length(params.time_start:.1:params.time_end), 1)';
T_SI_3_med = zeros(length(params.time_start:.1:params.time_end), 1)';
T_SI_3_high = zeros(length(params.time_start:.1:params.time_end), 1)';

for j = 1:length(params.time_start:.1:params.time_end)
   
   viremia_PI_1_low(j) =  quantile(viremia_1_PI(:,j), 0.0250); 
   viremia_PI_1_med(j) = quantile(viremia_1_PI(:,j), 0.5); 
   viremia_PI_1_high(j) = quantile(viremia_1_PI(:,j), 0.975); 
   
   viremia_SI_1_low(j) =  quantile(viremia_1_SI(:,j), 0.0250); 
   viremia_SI_1_med(j) = quantile(viremia_1_SI(:,j), 0.5); 
   viremia_SI_1_high(j) = quantile(viremia_1_SI(:,j), 0.975); 
   
   viremia_PI_2_low(j) =  quantile(viremia_2_PI(:,j), 0.0250); 
   viremia_PI_2_med(j) = quantile(viremia_2_PI(:,j), 0.5); 
   viremia_PI_2_high(j) = quantile(viremia_2_PI(:,j), 0.975); 
   
   viremia_SI_2_low(j) =  quantile(viremia_2_SI(:,j), 0.0250); 
   viremia_SI_2_med(j) = quantile(viremia_2_SI(:,j), 0.5); 
   viremia_SI_2_high(j) = quantile(viremia_2_SI(:,j), 0.975); 
   
   viremia_PI_3_low(j) =  quantile(viremia_3_PI(:,j), 0.0250); 
   viremia_PI_3_med(j) = quantile(viremia_3_PI(:,j), 0.5); 
   viremia_PI_3_high(j) = quantile(viremia_3_PI(:,j), 0.975); 
   
   viremia_SI_3_low(j) =  quantile(viremia_3_SI(:,j), 0.0250); 
   viremia_SI_3_med(j) = quantile(viremia_3_SI(:,j), 0.5); 
   viremia_SI_3_high(j) = quantile(viremia_3_SI(:,j), 0.975); 
   
   T_SI_1_low(j) =  quantile(T_1_SI(:,j), 0.0250); 
   T_SI_1_med(j) = quantile(T_1_SI(:,j), 0.5); 
   T_SI_1_high(j) = quantile(T_1_SI(:,j), 0.975); 
   
   T_SI_2_low(j) =  quantile(T_2_SI(:,j), 0.0250); 
   T_SI_2_med(j) = quantile(T_2_SI(:,j), 0.5); 
   T_SI_2_high(j) = quantile(T_2_SI(:,j), 0.975); 
   
   T_SI_3_low(j) =  quantile(T_3_SI(:,j), 0.0250); 
   T_SI_3_med(j) = quantile(T_3_SI(:,j), 0.5); 
   T_SI_3_high(j) = quantile(T_3_SI(:,j), 0.975); 
   
end

   subplot(x,y,1)
   shadedplot(TP1, viremia_PI_1_low,  viremia_PI_1_high, colors(1,:), 'none');
   hold on; 
   plot(TP1, viremia_PI_1_med, 'k', 'LineWidth', 2); 
   
   subplot(x,y,2)
   shadedplot(TP1, viremia_PI_2_low,  viremia_PI_2_high, colors(1,:), 'none');
   hold on; 
   plot(TP1, viremia_PI_2_med, 'k', 'LineWidth', 2); 
      
   subplot(x,y,3)
   shadedplot(TP1, viremia_PI_3_low,  viremia_PI_3_high, colors(1,:), 'none');
   hold on; 
   plot(TP1, viremia_PI_3_med, 'k', 'LineWidth', 2); 
      
   subplot(x,y,4)
   shadedplot(TP1, viremia_SI_1_low,  viremia_SI_1_high, colors(3,:), 'none');
   hold on; 
   plot(TP1, viremia_SI_1_med, 'k', 'LineWidth', 2); 
   
   subplot(x,y,5)
   shadedplot(TP1, viremia_SI_2_low,  viremia_SI_2_high, colors(3,:), 'none');
   hold on; 
   plot(TP1, viremia_SI_2_med, 'k', 'LineWidth', 2);    
   
   subplot(x,y,6)
   shadedplot(TP1, viremia_SI_3_low,  viremia_SI_3_high, colors(3,:), 'none');
   hold on; 
   plot(TP1, viremia_SI_3_med, 'k', 'LineWidth', 2); 
   
   fig.Units = 'inches';
   fig.PaperPosition = [1 1 7.5 5];
   fig.Position = fig.PaperPosition;

   print(fig,'fig6.tiff','-dtiff','-r300')

   
%plot T-cell dynamics
colors = colormap(cbrewer('qual', 'Paired', 10));

fig = figure();
x = 2; y = 3;

%model 3
load('chain_SS_3_IC_1_final') 
load('my_posterior_SS_3_IC_1_final')
my_chain = chain_multi_level;
my_posterior = my_posterior_multi_level; 

T2_1_SI = zeros(k, length(params.time_start:.1:params.time_end));
T2_2_SI = zeros(k, length(params.time_start:.1:params.time_end));
T2_3_SI = zeros(k, length(params.time_start:.1:params.time_end));

for i = 1:length(samples)
    params.beta = my_chain(samples(i), 1); 
    params.kappa  = my_chain(samples(i), 2); 
    params.q = my_chain(samples(i), 3); 
    params.Vinit = 10^(my_chain(samples(i), 8));
    params.sigma = my_chain(samples(i), 4);
    params.qT_1  = my_chain(samples(i), 5); 
    params.qT_2 = my_chain(samples(i), 5) + my_chain(samples(i), 6); 
    params.qT_3 = my_chain(samples(i), 5) + my_chain(samples(i), 7);
    
   params.qT = params.qT_1;
   [TS1, YS1] = ode45(@(t,y)SI(t, y, params), params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);        

   T2_1_SI(i,:) = YS3(:,5); 

   params.qT = params.qT_2; 
   [TS2, YS2] = ode45(@(t,y)SI(t, y, params), params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);        

   T2_2_SI(i,:) = YS3(:,5); 
   
   params.qT = params.qT_3; 
   [TS3, YS3] = ode45(@(t,y)SI(t, y, params), params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);        
          
   T2_3_SI(i,:) = YS3(:,5); 
   
end

T2_SI_1_low = zeros(length(params.time_start:.1:params.time_end), 1)';
T2_SI_1_med = zeros(length(params.time_start:.1:params.time_end), 1)';
T2_SI_1_high = zeros(length(params.time_start:.1:params.time_end), 1)';

T2_SI_2_low = zeros(length(params.time_start:.1:params.time_end), 1)';
T2_SI_2_med = zeros(length(params.time_start:.1:params.time_end), 1)';
T2_SI_2_high = zeros(length(params.time_start:.1:params.time_end), 1)';

T2_SI_3_low = zeros(length(params.time_start:.1:params.time_end), 1)';
T2_SI_3_med = zeros(length(params.time_start:.1:params.time_end), 1)';
T2_SI_3_high = zeros(length(params.time_start:.1:params.time_end), 1)';

for j = 1:length(params.time_start:.1:params.time_end)

   T2_SI_1_low(j) =  quantile(T2_1_SI(:,j), 0.0250); 
   T2_SI_1_med(j) = quantile(T2_1_SI(:,j), 0.5); 
   T2_SI_1_high(j) = quantile(T2_1_SI(:,j), 0.975); 
   
   T2_SI_2_low(j) =  quantile(T2_2_SI(:,j), 0.0250); 
   T2_SI_2_med(j) = quantile(T2_2_SI(:,j), 0.5); 
   T2_SI_2_high(j) = quantile(T2_2_SI(:,j), 0.975); 
   
   T2_SI_3_low(j) =  quantile(T2_3_SI(:,j), 0.0250); 
   T2_SI_3_med(j) = quantile(T2_3_SI(:,j), 0.5); 
   T2_SI_3_high(j) = quantile(T2_3_SI(:,j), 0.975); 
   
end

    subplot(x,y,1)
   shadedplot(TP1, T_SI_1_low, T_SI_1_high, colors(1,:), 'none');
   hold on; 
   plot(TP1', T_SI_1_med,'k', 'LineWidth', 2); 
   xlabel('time since infection (days)')
   ylabel('T-cells (cells/ml)');
   title('DENV-1')
  % set(gca, 'FontSize', 16);
   text(0.1,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
   xlim([0, 15]);   ylim([0, 2.2e7])

   subplot(x,y,2)
   shadedplot(TP1,T_SI_2_low, T_SI_2_high, colors(1,:), 'none');
   hold on; 
   plot(TP1', T_SI_2_med, 'k', 'LineWidth', 2); 
   xlabel('time since infection (days)')
   title('DENV-2')
   text(0.1,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top');
   xlim([0, 15]);   ylim([0, 2.2e7])

   subplot(x,y,3)
   shadedplot(TP1,T_SI_3_low, T_SI_3_high, colors(1,:), 'none');
   hold on; 
   plot(TP1', T_SI_3_med,'k', 'LineWidth', 2); 
   title('DENV-3')
   xlabel('time since infection (days)')
   text(0.1,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 
   xlim([0, 15]);   ylim([0, 2.2e7])

   subplot(x,y,4)
   shadedplot(TP1,T2_SI_1_low, T2_SI_1_high, colors(2,:), 'none');
   hold on; 
   plot(TP1', T2_SI_1_med,'k', 'LineWidth', 2); 
   xlabel('time since infection (days)')
   text(0.1,0.9,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 
   xlim([0, 15]);
   ylabel('T-cells (cells/ml)')
   ylim([0, 2.2e7])

   subplot(x,y,5)
   shadedplot(TP1,T2_SI_2_low, T2_SI_2_high, colors(2,:), 'none');
   hold on; 
   plot(TP1', T2_SI_2_med,'k', 'LineWidth', 2); 
   xlabel('time since infection (days)')
   text(0.1,0.9,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 
   xlim([0, 15]);
   ylim([0, 2.2e7])

   subplot(x,y,6)
   shadedplot(TP1,T2_SI_3_low, T2_SI_3_high, colors(2,:), 'none');
   hold on; 
   plot(TP1', T2_SI_3_med,'k', 'LineWidth', 2); 
   xlabel('time since infection (days)')
   text(0.1,0.9,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top');
   xlim([0, 15]);
   ylim([0, 2.2e7])

   fig.Units = 'inches';
   fig.PaperPosition = [1 1 7.5 5];
   fig.Position = fig.PaperPosition;

   print(fig,'figS6.tiff','-dtiff','-r300')

