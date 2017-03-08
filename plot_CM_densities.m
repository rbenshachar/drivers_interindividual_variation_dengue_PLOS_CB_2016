%plot_CM_densities

%FIGURE 4 IN PLOS CB paper
%plotting kernel densities of parameters of interest for several models

close all;

fig = figure; 

load('params');

colors = colormap(cbrewer('qual', 'Dark2', 5));
l = 150000; m = 100; n = 300000; 

    load('chain_OAS_IC_1_final')    
    my_chain = chain_CM; 
        
    x = 2; y = 2; 
    
    subplot(x, y,1)
    [f, xi] = ksdensity(params.deltaT + my_chain(l:m:n, 6));
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 2); %DF
    hold on; 
    line([params.deltaT, params.deltaT], [0, 6e6], 'Color', colors(3,:), 'LineWidth', 3);
    text(0.3,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
        ylim([0, 7e6])
    xlabel('\delta_T (day^{-1})'); ylabel('density')
    
   load('chain_OAS_alt_IC_1_final')
    my_chain = chain_CM; 
   
    hold on; 
   [f, xi] = ksdensity(params.deltaT + my_chain(l:m:n, 6)); 
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 2, 'LineStyle', '-.'); %DF 
    xlim([.5e-6, 4e-6])
    
     subplot(x, y, 2)
    [f, xi] = ksdensity(my_chain(l:m:n, 5));
    plot(xi, f, 'Color', colors(1,:), 'LineWidth', 2); %DF
    hold on; 
    [f, xi] = ksdensity(my_chain(l:m:n, 5) + my_chain(l:m:n, 7));
    plot(xi, f, 'Color', colors(3,:), 'LineWidth', 2); %DHF
    text(0.3,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 
    xlabel('q_T (day^{-1})'); ylabel('density');
    xlim([.5e-6, 4e-6])
    ylim([0, 7e6])
    legend('2^° DF', '2^° DHF')
    
    load('chain_ADE_IC_1_final')
    my_chain = chain_CM;

    subplot(x, y,3)
    [f, xi] = ksdensity(my_chain(l:m:n, 1));
    plot(xi, f, 'Color',colors(4,:), 'LineWidth', 2); 
    hold on;  
    [f, xi] = ksdensity(my_chain(l:m:n, 1) + my_chain(l:m:n, 6));
    text(0.1,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); 
    plot(xi, f, 'Color', colors(5,:),'LineWidth', 2); 
    xlabel('\beta ((copies/ml)^{-1} day^{-1})'); ylabel('density');
    ylim([0, 14e9])

  load('chain_ADE_IC_1_peak')   
  my_chain = chain_CM; 

    subplot(x, y, 4)
    [f, xi] = ksdensity(my_chain(l:n, 1));
    plot(xi, f, 'Color', colors(4,:), 'LineWidth', 2); 
    hold on;  
    [f, xi] = ksdensity(my_chain(l:n, 1) + my_chain(l:n, 7));
    plot(xi, f, 'Color', colors(5,:), 'LineWidth', 2); 
    xlabel('\beta ((copies/ml)^{-1} day^{-1})'); ylabel('density')
    legend('DF', '2^° DHF')
    ylim([0, 14e9])
    text(0.1,0.9,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top');

    fig.Units = 'inches';
    fig.PaperPosition = [1 1 7 6];
    fig.Position = fig.PaperPosition;

    print(fig,'fig4.tiff','-dtiff','-r400')
