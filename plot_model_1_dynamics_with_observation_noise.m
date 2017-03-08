%plot_model_1_dynamics_with_observation_noise

%FIGURE 2 in PLOS CB paper

%model 1

close all;

load('params'); %run save_params to get the parameters

colors = colormap(cbrewer('qual', 'Dark2', 5)); %using cbrewer for colors

load('chain_1_IC_1_final')
load('my_posterior_1_IC_1_final');

l = 150000; m = 100; n = 300000; 
my_chain = chain; 
params.time_end = 15; 
    
temp  = median((l:m:n));   

 params.beta = my_chain(temp(1), 1); 
 params.kappa = my_chain(temp(1), 2); 
 params.q = my_chain(temp(1), 3); 
 params.sigma = my_chain(temp(1), 4);
 params.qT = my_chain(temp(1), 5);
 params.Vinit = 10^(my_chain(temp(1), 6));
 params.IP = 5.9; 

 params.time_end = 15; 
     
  %PI
[T_best_PI, Y_best_PI] = ode45(@(t,y)PI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);  
[T_best_SI, Y_best_SI] = ode45(@(t,y)SI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);  

%IP_vals calculated from function sample_model_1_IP
load('IP_val_PI_1'); 
load('IP_val_SI_DF_1');
load('IP_val_SI_DHF_1');

%inter-individual variation
x = 3; y = 2; 

fig = figure;

subplot(x,y,1)
for i = 1:size(params.time_PI, 1)
    plot(IP_val_PI(i) + params.time_PI(i,:),  params.data_PI(i,:),  'k', 'LineWidth', 2);
    hold on; 
end
plot(0:15, 3.1761.*ones(16,1), 'k-.'); 
hold on; 
plot(0:15, 4.1761.*ones(16,1), 'k-.');
hold on; 
plot(T_best_PI, log10(Y_best_PI(:,3)),  'Color', colors(1,:), 'LineWidth', 2); 
xlim([0, 15]);ylim([0, 12.5])
text(0.1,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 
xlabel('time since infection (days)'); ylabel('viral load (log copies/ml)')

subplot(x,y,2)
for i = 1:size(params.time_SI_DF, 1)
    plot(IP_val_SI_DF(i) + params.time_SI_DF(i,:),  params.data_SI_DF(i,:),  'k', 'LineWidth', 2);
    hold on; 
end

for i = 1:size(params.time_SI_DHF, 1)
    plot(IP_val_SI_DHF(i) + params.time_SI_DHF(i,:),  params.data_SI_DHF(i,:),  'k', 'LineWidth', 2);
    hold on; 
end
plot(0:15, 3.1761.*ones(16,1), 'k-.'); 
hold on; 
plot(0:15, 4.1761.*ones(16,1), 'k-.');
hold on; 
plot(T_best_SI, log10(Y_best_SI(:,3)),  'Color', colors(3,:), 'LineWidth', 2); 
xlim([0, 15]); ylim([0, 12.5])
text(0.1,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
xlabel('time since infection (days)'); ylabel('viral load (log copies/ml)')

subplot(x, y, 3)
h = histogram(IP_val_PI); 
h.Normalization = 'probability';
h.NumBins = 20;
set(h,'facecolor', colors(1,:))
text(0.1,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
xlabel('incubation period (days)'); ylabel('density')

IP_val_SI = cat(1, IP_val_SI_DF, IP_val_SI_DHF); 

subplot(x,y,4) 
h1 = histogram(IP_val_SI_DF);
h1.Normalization = 'probability'; 
h1.NumBins = 20;
DF_vals = h1.Values./138; 
edges = h1.BinEdges(1:20);
h2 = histogram(IP_val_SI_DHF); 
h2.Normalization = 'probability'; 
h2.NumBins = 20;
h = bar(edges, [DF_vals; h2.Values./64]', 'Stacked');
h(1).FaceColor = colors(3,:); 
h(2).FaceColor = colors(4,:); 
legend('DF', 'DHF')
text(0.1,0.9,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %'FontSize', 16)
xlabel('incubation period (days)'); ylabel('density')

%individual trajectories with observation noise
samples = round(l + ((n-l))*rand(10, 1));
noise_PI = zeros(length(params.time_start:1:params.time_end), 1);
noise_SI = noise_PI; 

for i = 1:length(samples)
    params.beta = my_chain(samples(i), 1); 
    params.kappa = my_chain(samples(i), 2); 
    params.q = my_chain(samples(i), 3); 
    params.sigma = my_chain(samples(i), 4); 
    params.qT = my_chain(samples(i), 5); 
    params.Vinit = 10^(my_chain(samples(i), 6)); 
    
     [T, Y] = ode45(@(t,y)PI(t, y, params),params.time_start:1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);  
     
     V_PI = log10(Y(:,3)); 
     for j = 1:length(noise_PI)
         noise_PI(j) = normrnd(0,1); 
     end
     
     V_with_noise_PI = V_PI + noise_PI; 
     
     subplot(x,y,5)
     plot(T, V_with_noise_PI,  'Color', colors(1,:))
     hold on; 
     
     [T2, Y2] = ode45(@(t,y)SI(t, y, params),params.time_start:1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);   
     V_SI = log10(Y2(:,3));
     
     for j = 1:length(noise_SI)
         noise_SI(j) = normrnd(0,1); 
     end
     V_with_noise_SI = V_SI + noise_SI; 
     
     subplot(x,y,6)
     plot(T, V_with_noise_SI,  'Color', colors(3,:))
     hold on; 

end

subplot(x,y,5)
text(0.1,0.9,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 
ylabel('viral load (log copies/ml)')
xlabel('time since infection (days)')
xlim([0, 15]); ylim([0, 12.5]);

subplot(x,y,6)
text(0.1,0.9,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top');
ylabel('viral load (log copies/ml)')
xlabel('time since infection (days)')
xlim([0, 15]); ylim([0, 12.5])

 fig.Units = 'inches';
    fig.PaperPosition = [1 1 6 6];
    fig.Position = fig.PaperPosition;

print(fig,'fig2.tiff','-dtiff','-r300')

