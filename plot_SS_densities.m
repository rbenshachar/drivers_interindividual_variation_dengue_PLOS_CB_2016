%plot_SS_densities

%FIGURE 5 IN PLOS CB paper
%plotting kernel densities of parameters of interest for several models

close all;

fig = figure; 

load('params');
l = 150000; m = 100; n = 300000; 
colors = colormap(cbrewer('qual', 'Paired', 10));

load('chain_SS_1_IC_1_final')
my_chain = chain_multi_level; 
    
    x = 2; y = 2; 
    subplot(x, y,1)
    [f, xi] = ksdensity(my_chain(l:m:n, 1));
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 3); %DENV-1
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 6));
    plot(xi, f, 'Color', colors(2,:), 'LineWidth', 3); %DENV-2
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 7));
    plot(xi, f, 'Color', colors(3,:), 'LineWidth', 3); %DENV-3
    xlabel('\beta ((copies/ml)^{-1} day^{-1})'); ylabel('density');

    load('chain_SS_2_IC_1_final')    
    my_chain = chain_multi_level; 
    
    subplot(x, y,2)
    [f, xi] = ksdensity(my_chain(l:m:n, 3));
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 3); %DENV-1
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 3) + my_chain(l:m:n, 6));
    plot(xi, f, 'Color', colors(2,:), 'LineWidth',3); %DENV-2
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 3) + my_chain(l:m:n, 7));
    plot(xi, f, 'Color', colors(3,:), 'LineWidth', 3); %DENV-3
    xlabel('q (day^{-1})'); ylabel('density')
    text(0.7,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 

    load('chain_SS_3_IC_1_final')    
    my_chain = chain_multi_level;
    
    subplot(x, y,3)
    [f, xi] = ksdensity(my_chain(l:m:n, 5));
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 3); %DENV-1
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 5) + my_chain(l:m:n, 6));
    plot(xi, f, 'Color', colors(2,:), 'LineWidth', 3); %DENV-2
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 5) + my_chain(l:m:n, 7));
    plot(xi, f, 'Color', colors(3,:), 'LineWidth', 3); %DENV-3
    xlabel('q_T ((cells/ml)^{-1} day^{-1})'); ylabel('density')
    text(0.7,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top');
    
    load('chain_SS_1_IC_1_peak_final_prior'); 
    my_chain = chain_multi_level; 
    
    subplot(x, y,4)
    [f, xi] = ksdensity(my_chain(l:m:n, 1));
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 3); %DENV-1
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 7));
    plot(xi, f, 'Color', colors(2,:), 'LineWidth', 3); %DENV-2
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 8));
    plot(xi, f, 'Color', colors(3,:), 'LineWidth', 3); %DENV-3
    xlabel('\beta ((copies/ml)^{-1} day^{-1})'); ylabel('density')
    text(0.7,0.9,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 
    legend('DENV-1', 'DENV-2', 'DENV-3'); 

fig.Units = 'inches';
fig.PaperPosition = [1 1 7 6];
fig.Position = fig.PaperPosition;

print(fig,'fig5.tiff','-dtiff','-r300')


