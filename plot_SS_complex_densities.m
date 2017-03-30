%plot_SS_complex_densities

%FIGURE S7

close all;

load('params');

fig = figure();

parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20, 1e-3, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0, 1];

l = 150000; m = 1; n = 300000; 
colors = colormap(cbrewer('qual', 'Paired', 10));

x = 1; y = 2; 
 
    load('chain_SS_1_ADE_IC_1_final_complex')   
    my_chain = chain_multi_level; 
    
    subplot(x, y,1)
    [f, xi] = ksdensity(my_chain(l:m:n, 1));
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 3); %DF
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 6));
    plot(xi, f, 'Color', colors(2,:), 'LineWidth', 3); %DF
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 7));
    plot(xi, f, 'Color', colors(3,:), 'LineWidth', 3); %D
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 8));
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 3, 'LineStyle', '-.'); %DF
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 6) + my_chain(l:m:n, 9));
    plot(xi, f, 'Color', colors(2,:), 'LineWidth', 3, 'LineStyle', '-.'); %DF
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 7) + my_chain(l:m:n, 10));
    plot(xi, f, 'Color', colors(3,:), 'LineWidth', 3, 'LineStyle', '-.'); 
    xlabel('\beta ((copies/ml)^{-1} day^{-1})'); ylabel('density')
    text(0.7,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
    xlim([2e-10, 10e-10])
    %set(gca, 'FontSize', 16)
        
    load('chain_SS_1_ADE_IC_1_peak_final_prior_complex')   
    my_chain = chain_multi_level; 
    
    subplot(x, y,2)
    [f, xi] = ksdensity(my_chain(l:m:n, 1));
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 3); %DF
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 7));
    plot(xi, f, 'Color', colors(2,:), 'LineWidth', 3); %DF
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 8));
    plot(xi, f, 'Color', colors(3,:), 'LineWidth', 3); %D
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 9));
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 3, 'LineStyle', '-.'); %DF
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 7) + my_chain(l:m:n, 10));
    plot(xi, f, 'Color', colors(2,:), 'LineWidth', 3, 'LineStyle', '-.'); %DF
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1)  + my_chain(l:m:n, 8) + my_chain(l:m:n, 11));
    plot(xi, f, 'Color', colors(3,:), 'LineWidth', 3, 'LineStyle', '-.'); 
    xlabel('\beta ((copies/ml)^{-1} day^{-1})'); ylabel('density')
    text(0.7,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
    xlim([2e-10, 10e-10])
    %set(gca, 'FontSize', 16)
       

   fig.Units = 'inches';
   fig.PaperPosition = [1 1 7.5 4];
   fig.Position = fig.PaperPosition;

   print(fig,'figS8.tiff','-dtiff','-r300')




