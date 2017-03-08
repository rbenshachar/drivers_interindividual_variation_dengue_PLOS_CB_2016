%plot_model_1_dynamics_with_variance

%FINAL FIGURE 3 in PLOS CB paper
%model 1 dynamics including variance

close all;
clear all; 

fig = figure; 

load('params');
colors = colormap(cbrewer('qual', 'Dark2', 3));

parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20, 1e-3, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0, 1];

load('chain_1_IC_1_final')
load('my_posterior_1_IC_1_final');

l = 150000; m = 100; n = 300000; 
my_chain = chain; 

k = 100; 
samples = round(l + ((n-l))*rand(k, 1));

params.time_end = 15; 
params.IP = 5.9; 
     
uninfected_cells_PI = zeros(k, length(params.time_start:.01:params.time_end));
NK_cells_PI = zeros(k, length(params.time_start:.01:params.time_end)); 
viremia_PI = zeros(k, length(params.time_start:.01:params.time_end)); 
uninfected_cells_SI = zeros(k, length(params.time_start:.01:params.time_end));
viremia_SI = zeros(k, length(params.time_start:.01:params.time_end));
NK_cells_SI = zeros(k, length(params.time_start:.01:params.time_end));
T_cells_SI = zeros(k, length(params.time_start:.01:params.time_end)); 

 x = 2; y = 2;

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
     uninfected_cells_PI(i,:) = log10(Y(:,1)); 
     viremia_PI(i,:) = log10(Y(:,3)); 
     NK_cells_PI(i,:) = Y(:,4); 
     
     T_SI = T2; 
     uninfected_cells_SI(i,:) = log10(Y2(:,1)); 
     viremia_SI(i,:) = log10(Y2(:,3)); 
     NK_cells_SI(i,:) = Y2(:,4); 
     T_cells_SI(i,:) = Y2(:,5); 

end

uninfected_cells_PI_low = zeros(length(params.time_start:.01:params.time_end), 1)';
uninfected_cells_PI_med = zeros(length(params.time_start:.01:params.time_end), 1)';
uninfected_cells_PI_high = zeros(length(params.time_start:.01:params.time_end), 1)';

uninfected_cells_SI_low = zeros(length(params.time_start:.01:params.time_end), 1)';
uninfected_cells_SI_med = zeros(length(params.time_start:.01:params.time_end), 1)';
uninfected_cells_SI_high = zeros(length(params.time_start:.01:params.time_end), 1)';

viremia_PI_low = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_PI_med = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_PI_high = zeros(length(params.time_start:.01:params.time_end), 1)';

viremia_SI_low = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_SI_med = zeros(length(params.time_start:.01:params.time_end), 1)';
viremia_SI_high = zeros(length(params.time_start:.01:params.time_end), 1)';

NK_cells_PI_low = zeros(length(params.time_start:.01:params.time_end), 1)';
NK_cells_PI_med = zeros(length(params.time_start:.01:params.time_end), 1)';
NK_cells_PI_high = zeros(length(params.time_start:.01:params.time_end), 1)';

NK_cells_SI_low = zeros(length(params.time_start:.01:params.time_end), 1)';
NK_cells_SI_med = zeros(length(params.time_start:.01:params.time_end), 1)';
NK_cells_SI_high = zeros(length(params.time_start:.01:params.time_end), 1)';

T_cells_SI_low = zeros(length(params.time_start:.01:params.time_end), 1)';
T_cells_SI_med = zeros(length(params.time_start:.01:params.time_end), 1)';
T_cells_SI_high = zeros(length(params.time_start:.01:params.time_end), 1)';

for j = 1:length(params.time_start:.01:params.time_end)
   uninfected_cells_PI_low(j) =  quantile(uninfected_cells_PI(:,j), 0.0250); 
   uninfected_cells_PI_med(j) =  quantile(uninfected_cells_PI(:,j), 0.50); 
   uninfected_cells_PI_high(j) = quantile(uninfected_cells_PI(:,j), 0.975); 
   
   uninfected_cells_SI_low(j) =  quantile(uninfected_cells_SI(:,j), 0.0250); 
   uninfected_cells_SI_med(j) =  quantile(uninfected_cells_SI(:,j), 0.50); 
   uninfected_cells_SI_high(j) = quantile(uninfected_cells_SI(:,j), 0.975); 
   
   viremia_PI_low(j) =  quantile(viremia_PI(:,j), 0.0250); 
   viremia_PI_med(j) = quantile(viremia_PI(:,j), 0.5); 
   viremia_PI_high(j) = quantile(viremia_PI(:,j), 0.975); 
   
   viremia_SI_low(j) =  quantile(viremia_SI(:,j), 0.0250); 
   viremia_SI_med(j) = quantile(viremia_SI(:,j), 0.5); 
   viremia_SI_high(j) = quantile(viremia_SI(:,j), 0.975); 
   
   NK_cells_PI_low(j) = quantile(NK_cells_PI(:,j), .0250); 
   NK_cells_PI_med(j) = quantile(NK_cells_PI(:,j), .50); 
   NK_cells_PI_high(j) = quantile(NK_cells_PI(:,j), .9750); 
   
   NK_cells_SI_low(j) = quantile(NK_cells_SI(:,j), .0250); 
   NK_cells_SI_med(j) = quantile(NK_cells_SI(:,j), .50); 
   NK_cells_SI_high(j) = quantile(NK_cells_SI(:,j), .9750); 
   
   T_cells_SI_low(j) = quantile(T_cells_SI(:,j), .0250); 
   T_cells_SI_med(j) = quantile(T_cells_SI(:,j), 0.5); 
   T_cells_SI_high(j) = quantile(T_cells_SI(:,j), .9750); 
   
end

subplot(x,y,1) 
shadedplot(T_PI, uninfected_cells_PI_low,  uninfected_cells_PI_high, colors(1,:), 'none');
hold on;
shadedplot(T_SI, uninfected_cells_SI_low,  uninfected_cells_SI_high, colors(3,:), 'none');
hold on; 
h3 = plot(T_PI, uninfected_cells_PI_med, 'k', 'LineWidth', 2); 
hold on; 
h4 = plot(T_SI, uninfected_cells_SI_med, 'k-.', 'LineWidth', 2);
text(0.1,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top');
legend([h3, h4], {'1^°', '2^°'})
xlabel('time since infection (days)')
ylabel('uninfected cells (log cells/ml)') 
ylim([3,7])

subplot(x,y,2) 
shadedplot(T_PI, viremia_PI_low,  viremia_PI_high, colors(1,:), 'none');
hold on;
shadedplot(T_SI, viremia_SI_low,  viremia_SI_high, colors(3,:), 'none');
hold on; 
plot(T_PI, viremia_PI_med, 'k', 'LineWidth', 2); 
hold on; 
plot(T_SI, viremia_SI_med,  'k-.', 'LineWidth', 2);
text(0.1,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 
xlabel('time since infection (days)')
ylabel('viral load (log copies/ml)')
ylim([-4, 12])

subplot(x,y,3) 
shadedplot(T_PI, NK_cells_PI_low,  NK_cells_PI_high, colors(1,:), 'none');
hold on;
shadedplot(T_SI, NK_cells_SI_low,  NK_cells_SI_high, colors(3,:), 'none');
hold on; 
plot(T_PI, NK_cells_PI_med,  'k', 'LineWidth', 2); 
hold on; 
plot(T_SI, NK_cells_SI_med,  'k-.', 'LineWidth', 2);
text(0.1,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 
xlabel('time since infection (days)')
ylabel('NK cells (cells/ml)')
ylim([-5, 4000])
    
subplot(x,y,4)  
shadedplot(T_SI, T_cells_SI_low,  T_cells_SI_high, colors(3,:), 'none');
hold on; 
plot(T_SI, T_cells_SI_med,  'k-.',  'LineWidth', 2);
text(0.1,0.9,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top');
xlabel('time since infection (days)')
ylabel('T-cells (cells/ml)')
ylim([-5,5e6])

fig.Units = 'inches';
fig.PaperPosition = [1 1 6 6];
fig.Position = fig.PaperPosition;

print(fig,'fig3.tiff','-dtiff','-r300')
