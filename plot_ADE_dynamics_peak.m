%plot_ADE_dynamics_peak

%FIGURE S5 IN PLOS CB paper
%plots model dynamics of model ADE peak, stratifying by primary DF,
%secondary DF and secondary DHF

close all;
load('params')

fig = figure();

colors = colormap(cbrewer('qual', 'Set2', 3));

data_with_peak; %identifying subset of individuals with viral peak and puts this data in params structure

load('chain_ADE_IC_1_peak')
load('my_posterior_ADE_IC_1_peak')

my_posterior = my_posterior_CM;

l = 150000; m = 1; n = 300000; 
my_chain = chain_CM; 

k = 100; 
samples = round(l + ((n-l))*rand(k, 1));

params.time_end = 15; 

params.IP = 5.9; 
params.Vinit = 10^(-3.2); 
     
uninfected_cells_PI = zeros(k, length(params.time_start:.1:params.time_end));
NK_cells_PI = zeros(k, length(params.time_start:.1:params.time_end)); 
viremia_PI = zeros(k, length(params.time_start:.1:params.time_end)); 

uninfected_cells_SI_DF = zeros(k, length(params.time_start:.1:params.time_end));
viremia_SI_DF = zeros(k, length(params.time_start:.1:params.time_end));
NK_cells_SI_DF = zeros(k, length(params.time_start:.1:params.time_end));
T_cells_SI_DF = zeros(k, length(params.time_start:.1:params.time_end)); 

uninfected_cells_SI_DHF = zeros(k, length(params.time_start:.1:params.time_end));
viremia_SI_DHF = zeros(k, length(params.time_start:.1:params.time_end));
NK_cells_SI_DHF = zeros(k, length(params.time_start:.1:params.time_end));
T_cells_SI_DHF = zeros(k, length(params.time_start:.1:params.time_end)); 

x = 2; y = 2;

for i = 1:length(samples)
    
    params.beta = my_chain(samples(i), 1); 
    params.kappa = my_chain(samples(i), 2); 
    params.q = my_chain(samples(i), 3); 
    params.IP = my_chain(samples(i), 4); 
    params.sigma = my_chain(samples(i), 5); 
    params.qT = my_chain(samples(i), 6); 
        
    [T, Y] = ode45(@(t,y)PI(t, y, params),params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);  
    [T2, Y2] = ode45(@(t,y)SI(t, y, params),params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);     
    
    params.beta = my_chain(samples(i), 1) + my_chain(samples(i), 7); 
    [T3, Y3] = ode45(@(t,y)SI(t, y, params),params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);     
    
    uninfected_cells_PI(i,:) = log10(Y(:,1)); 
    viremia_PI(i,:) = log10(Y(:,3)); 
    NK_cells_PI(i,:) = Y(:,4); 
    
    uninfected_cells_SI_DF(i,:) = log10(Y2(:,1)); 
    viremia_SI_DF(i,:) = log10(Y2(:,3)); 
    NK_cells_SI_DF(i,:) = Y2(:,4); 
    T_cells_SI_DF(i,:) = Y2(:,5); 
    
    uninfected_cells_SI_DHF(i,:) = log10(Y3(:,1)); 
    viremia_SI_DHF(i,:) = log10(Y3(:,3)); 
    NK_cells_SI_DHF(i,:) = Y3(:,4); 
    T_cells_SI_DHF(i,:) = Y3(:,5); 
    
end

uninfected_cells_PI_low = zeros(length(params.time_start:.1:params.time_end), 1)';
uninfected_cells_PI_med = zeros(length(params.time_start:.1:params.time_end), 1)';
uninfected_cells_PI_high = zeros(length(params.time_start:.1:params.time_end), 1)';

uninfected_cells_SI_DF_low = zeros(length(params.time_start:.1:params.time_end), 1)';
uninfected_cells_SI_DF_med = zeros(length(params.time_start:.1:params.time_end), 1)';
uninfected_cells_SI_DF_high = zeros(length(params.time_start:.1:params.time_end), 1)';

uninfected_cells_SI_DHF_low = zeros(length(params.time_start:.1:params.time_end), 1)';
uninfected_cells_SI_DHF_med = zeros(length(params.time_start:.1:params.time_end), 1)';
uninfected_cells_SI_DHF_high = zeros(length(params.time_start:.1:params.time_end), 1)';

viremia_PI_low = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_PI_med = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_PI_high = zeros(length(params.time_start:.1:params.time_end), 1)';

viremia_SI_DF_low = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_SI_DF_med = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_SI_DF_high = zeros(length(params.time_start:.1:params.time_end), 1)';

viremia_SI_DHF_low = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_SI_DHF_med = zeros(length(params.time_start:.1:params.time_end), 1)';
viremia_SI_DHF_high = zeros(length(params.time_start:.1:params.time_end), 1)';

NK_cells_PI_low = zeros(length(params.time_start:.1:params.time_end), 1)';
NK_cells_PI_med = zeros(length(params.time_start:.1:params.time_end), 1)';
NK_cells_PI_high = zeros(length(params.time_start:.1:params.time_end), 1)';

NK_cells_SI_DF_low = zeros(length(params.time_start:.1:params.time_end), 1)';
NK_cells_SI_DF_med = zeros(length(params.time_start:.1:params.time_end), 1)';
NK_cells_SI_DF_high = zeros(length(params.time_start:.1:params.time_end), 1)';

T_cells_SI_DF_low = zeros(length(params.time_start:.1:params.time_end), 1)';
T_cells_SI_DF_med = zeros(length(params.time_start:.1:params.time_end), 1)';
T_cells_SI_DF_high = zeros(length(params.time_start:.1:params.time_end), 1)';

NK_cells_SI_DHF_low = zeros(length(params.time_start:.1:params.time_end), 1)';
NK_cells_SI_DHF_med = zeros(length(params.time_start:.1:params.time_end), 1)';
NK_cells_SI_DHF_high = zeros(length(params.time_start:.1:params.time_end), 1)';

T_cells_SI_DHF_low = zeros(length(params.time_start:.1:params.time_end), 1)';
T_cells_SI_DHF_med = zeros(length(params.time_start:.1:params.time_end), 1)';
T_cells_SI_DHF_high = zeros(length(params.time_start:.1:params.time_end), 1)';

for j = 1:length(params.time_start:.1:params.time_end)
   uninfected_cells_PI_low(j) =  quantile(uninfected_cells_PI(:,j), 0.0250); 
   uninfected_cells_PI_med(j) =  quantile(uninfected_cells_PI(:,j), 0.50); 
   uninfected_cells_PI_high(j) = quantile(uninfected_cells_PI(:,j), 0.975); 
   
   uninfected_cells_SI_DF_low(j) =  quantile(uninfected_cells_SI_DF(:,j), 0.0250); 
   uninfected_cells_SI_DF_med(j) =  quantile(uninfected_cells_SI_DF(:,j), 0.50); 
   uninfected_cells_SI_DF_high(j) = quantile(uninfected_cells_SI_DF(:,j), 0.975); 
   
   uninfected_cells_SI_DHF_low(j) =  quantile(uninfected_cells_SI_DHF(:,j), 0.0250); 
   uninfected_cells_SI_DHF_med(j) =  quantile(uninfected_cells_SI_DHF(:,j), 0.50); 
   uninfected_cells_SI_DHF_high(j) = quantile(uninfected_cells_SI_DHF(:,j), 0.975); 
   
   viremia_PI_low(j) =  quantile(viremia_PI(:,j), 0.0250); 
   viremia_PI_med(j) = quantile(viremia_PI(:,j), 0.5); 
   viremia_PI_high(j) = quantile(viremia_PI(:,j), 0.975); 
   
   viremia_SI_DF_low(j) =  quantile(viremia_SI_DF(:,j), 0.0250); 
   viremia_SI_DF_med(j) = quantile(viremia_SI_DF(:,j), 0.5); 
   viremia_SI_DF_high(j) = quantile(viremia_SI_DF(:,j), 0.975); 
   
   viremia_SI_DHF_low(j) =  quantile(viremia_SI_DHF(:,j), 0.0250); 
   viremia_SI_DHF_med(j) = quantile(viremia_SI_DHF(:,j), 0.5); 
   viremia_SI_DHF_high(j) = quantile(viremia_SI_DHF(:,j), 0.975); 
   
   NK_cells_PI_low(j) = quantile(NK_cells_PI(:,j), .0250); 
   NK_cells_PI_med(j) = quantile(NK_cells_PI(:,j), .50); 
   NK_cells_PI_high(j) = quantile(NK_cells_PI(:,j), .9750); 
   
   NK_cells_SI_DF_low(j) = quantile(NK_cells_SI_DF(:,j), .0250); 
   NK_cells_SI_DF_med(j) = quantile(NK_cells_SI_DF(:,j), .50); 
   NK_cells_SI_DF_high(j) = quantile(NK_cells_SI_DF(:,j), .9750); 
   
   NK_cells_SI_DHF_low(j) = quantile(NK_cells_SI_DHF(:,j), .0250); 
   NK_cells_SI_DHF_med(j) = quantile(NK_cells_SI_DHF(:,j), .50); 
   NK_cells_SI_DHF_high(j) = quantile(NK_cells_SI_DHF(:,j), .9750); 
   
   T_cells_SI_DF_low(j) = quantile(T_cells_SI_DF(:,j), .0250); 
   T_cells_SI_DF_med(j) = quantile(T_cells_SI_DF(:,j), 0.5); 
   T_cells_SI_DF_high(j) = quantile(T_cells_SI_DF(:,j), .9750); 
   
   T_cells_SI_DHF_low(j) = quantile(T_cells_SI_DHF(:,j), .0250); 
   T_cells_SI_DHF_med(j) = quantile(T_cells_SI_DHF(:,j), 0.5); 
   T_cells_SI_DHF_high(j) = quantile(T_cells_SI_DHF(:,j), .9750); 
end

subplot(x,y,1)
shadedplot(T, uninfected_cells_PI_low,  uninfected_cells_PI_high, colors(1,:), 'none');
hold on;
shadedplot(T, uninfected_cells_SI_DF_low,  uninfected_cells_SI_DF_high, colors(2,:), 'none');
hold on; 
shadedplot(T, uninfected_cells_SI_DHF_low,  uninfected_cells_SI_DHF_high, colors(3,:), 'none');
hold on; 
plot(T, uninfected_cells_PI_med, 'k', 'LineWidth', 2); 
hold on; 
plot(T, uninfected_cells_SI_DF_med, 'k-.', 'LineWidth', 2);
hold on; 
plot(T, uninfected_cells_SI_DHF_med, 'k--', 'LineWidth', 2);
text(0.1,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
xlabel('time since infection (days)'); ylabel('uninfected cells (log cells/ml)')
%set(gca, 'FontSize', 16); 
ylim([3,7])

subplot(x,y,2)
shadedplot(T, viremia_PI_low, viremia_PI_high, colors(1,:), 'none');
hold on;
shadedplot(T, viremia_SI_DF_low, viremia_SI_DF_high, colors(2,:), 'none');
hold on; 
shadedplot(T, viremia_SI_DHF_low, viremia_SI_DHF_high, colors(3,:), 'none');
hold on; 
h4 = plot(T, viremia_PI_med, 'k', 'LineWidth', 2); 
hold on; 
h5 = plot(T, viremia_SI_DF_med, 'k-.', 'LineWidth', 2);
hold on; 
h6 = plot(T, viremia_SI_DHF_med, 'k--', 'LineWidth', 2);
text(0.1,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
%legend([h4 h5 h6], {'1^° DF', '2^° DF', '2^° DHF'})
xlabel('time since infection (days)'); ylabel('viral load (log copies/ml)')
%set(gca, 'FontSize', 16); 
ylim([-4,12])

subplot(x,y,3)
shadedplot(T, NK_cells_PI_low, NK_cells_PI_high, colors(1,:), 'none');
hold on;
shadedplot(T, NK_cells_SI_DF_low, NK_cells_SI_DF_high, colors(2,:), 'none');
hold on; 
shadedplot(T, NK_cells_SI_DHF_low, NK_cells_SI_DHF_high, colors(3,:), 'none');
hold on; 
plot(T, NK_cells_PI_med, 'k', 'LineWidth', 2); 
hold on; 
plot(T, NK_cells_SI_DF_med, 'k-.', 'LineWidth', 2);
hold on; 
plot(T, NK_cells_SI_DHF_med, 'k--', 'LineWidth', 2);
text(0.1,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
xlabel('time since infection (days)'); ylabel('NK cells (cells/ml)')
%set(gca, 'FontSize', 16); 
ylim([0, 5000])

subplot(x,y,4)
shadedplot(T, T_cells_SI_DF_low, T_cells_SI_DF_high, colors(2,:), 'none');
hold on; 
shadedplot(T, T_cells_SI_DHF_low, T_cells_SI_DHF_high, colors(3,:), 'none');
hold on; 
plot(T, T_cells_SI_DF_med, 'k-.', 'LineWidth', 2);
hold on; 
plot(T, T_cells_SI_DHF_med, 'k--', 'LineWidth', 2);
text(0.1,0.9,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
xlabel('time since infection (days)'); ylabel('T cells (cells/ml)')
%set(gca, 'FontSize', 16); 
ylim([0, 6e6])


   fig.Units = 'inches';
   fig.PaperPosition = [1 1 7.5 5];
   fig.Position = fig.PaperPosition;

   print(fig,'figS5.tiff','-dtiff','-r300')