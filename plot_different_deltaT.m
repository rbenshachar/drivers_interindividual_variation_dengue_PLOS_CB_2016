%plot_different_deltaT

%FIGURE S3 IN PLOS CB 2016 paper
%Plots viral load and T-cell dynamics when model 1 is fit for different
%values of deltaT

close all; 

fig = figure;

load('params');

params.IP = 5.9; 

params.time_end = 15;

    x = 1; y = 2;

    load('chain_1_IC_1_deltaT_5')
    load('my_posterior_1_IC_1_deltaT_5')
    params.deltaT = 1e-5;
    l = 150000; n = 300000; 
    my_chain = chain; 
    posterior_2 = median(my_posterior(l:n)); 
    temp = find(my_posterior == posterior_2); 
    
	params.beta  = my_chain(temp(1), 1) ;
    params.kappa = my_chain(temp(1), 2) ;
    params.q  = my_chain(temp(1), 3) ;
    params.Vinit =  10^my_chain(temp(1), 6) ;
    params.sigma = my_chain(temp(1), 4);
    params.qT = my_chain(temp(1), 5);
    
    
    [T, Y] = ode45(@(t,y)SI(t, y, params),params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit, params.Tinit]);  
    
    subplot(x, y, 2)
    semilogy(T, Y(:,5), 'm', 'LineWidth', 2); 
    xlabel('time since infection (days)')
    ylabel('T-cells (cells/ml)')
    text(0.1,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16);
    %set(gca, 'FontSize', 16); 
    hold on; 
   
    subplot(x, y, 1)
    plot(T, log10(Y(:,3)), 'm', 'LineWidth', 2); 
    xlabel('time since infection (days)')
    ylabel('viral load (log copies/ml)')
    text(0.1,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16);
    %set(gca, 'FontSize', 16); 
    hold on; 
    
    load('chain_1_IC_1_final')
    load('my_posterior_1_IC_1_final')
    params.deltaT = 1e-6;
    l = 150000; m = 100; n = 300000; 
    my_chain = chain; 
     posterior_3 = median(my_posterior(l:n)); 
    temp = find(my_posterior == posterior_3); 

    params.beta  = my_chain(temp(1), 1); 
    params.kappa = my_chain(temp(1), 2); 
    params.q  = my_chain(temp(1), 3); 
    params.Vinit =  10^my_chain(temp(1), 6); 
    params.sigma = my_chain(temp(1), 4); 
    params.qT = my_chain(temp(1), 5);    
    
    [T, Y] = ode45(@(t,y)SI(t, y, params),params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit, params.Tinit]);  
      
    subplot(x, y, 2)
    semilogy(T, Y(:,5), 'g', 'LineWidth', 2);
    hold on; 
    
    subplot(x, y,1)
    plot(T, log10(Y(:,3)), 'g', 'LineWidth', 2); 
    hold on; 

    load('chain_1_IC_1_deltaT_7')
    load('my_posterior_1_IC_1_deltaT_7')
    params.deltaT = 1e-7;
    l = 150000; m = 100; n = 300000; 
    my_chain = chain; 
    posterior_4 = median(my_posterior(l:n)); 
    temp = find(my_posterior == posterior_4); 

	params.beta  = my_chain(temp(1), 1); 
    params.kappa = my_chain(temp(1), 2); 
    params.q  = my_chain(temp(1), 3); 
    params.Vinit =  10^my_chain(temp(1), 6); 
    params.sigma = my_chain(temp(1), 4); 
    params.qT = my_chain(temp(1), 5);
    
    
    [T, Y] = ode45(@(t,y)SI(t, y, params),params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit, params.Tinit]);  

    subplot(x, y, 2)
    semilogy(T, Y(:,5), 'r', 'LineWidth', 2);
    hold on; 
     
    subplot(x, y,1)
    plot(T, log10(Y(:,3)), 'r', 'LineWidth', 2); 
    hold on;   
        
    load('chain_1_IC_1_deltaT_8')
    load('my_posterior_1_IC_1_deltaT_8')
    params.deltaT = 1e-8;
    l = 150000; m = 100; n = 300000; 
    my_chain = chain; 
    posterior_5 = median(my_posterior(l:n)); 
    temp = find(my_posterior == posterior_5); 
    
	params.beta  = my_chain(temp(1), 1); 
    params.kappa = my_chain(temp(1), 2); 
    params.q  = my_chain(temp(1), 3); 
    params.Vinit =  10^my_chain(temp(1), 6); 
    params.sigma = my_chain(temp(1), 4); 
    params.qT = my_chain(temp(1), 5);
    

    [T, Y] = ode45(@(t,y)SI(t, y, params),params.time_start:.1:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit, params.Tinit]);  

    subplot(x, y,1)
    plot(T, log10(Y(:,3)), 'c', 'LineWidth', 2); 
    hold on; 
    plot(0:15,3.1761.*ones(16,1), 'k-.'); 
    hold on; 
    plot(0:15,4.1761.*ones(16,1), 'k-.');
    legend('\delta_T = 10^{-5}/day', '10^{-6}/day', '10^{-7}/day', '10^{-8}/day', 'Location', 'south')
   
    subplot(x, y, 2)
    semilogy(T, Y(:,5), 'c', 'LineWidth', 2);
    
   fig.Units = 'inches';
   fig.PaperPosition = [1 1 7.5 3.5];
   fig.Position = fig.PaperPosition;

   print(fig,'figS3.tiff','-dtiff','-r300')



 