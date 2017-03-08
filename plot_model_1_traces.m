function void = plot_model_1_traces()

%SUPPLEMENTARY FIGURE 1 in PLOS CB 2016 paper
%plots traces of model 1

close all; 

fig = figure;
z = 6; 
names = {'\beta', '\kappa', 'q', '\sigma_I', 'q_T', 'log(V_0)','log-likelihood'};
caption_label = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)'};

l = 1; m = 1000; n = 300000; 

%chain with IC1
load('chain_1_IC_1_final')
load('my_posterior_1_IC_1_final')

for i = 1:z
 subplot(4,2, i)
    plot(l:m:n, chain(l:m:n, i), 'r'); 
    xlabel({'MCMC', 'iteration'}); ylabel(names(i));
    text(0.9,0.9,caption_label(i),'Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
    if i == 6
       ylim([-6, 0])
    end
    hold on;
end

subplot(4,2,z+1)
    plot(l:m:n, my_posterior(l:m:n), 'r'); 
    xlabel({'MCMC', 'iteration'}); ylabel('log-likelihood'); 
        text(0.9,0.9,caption_label(end),'Units', 'Normalized', 'VerticalAlignment', 'Top'); %, 'FontSize', 16)
        hold on;
  
%chain with IC2
load('chain_1_IC_2_final')
load('my_posterior_1_IC_2_final')
    
for i = 1:z
    subplot(4,2, i)
    plot(l:m:n, chain(l:m:n, i), 'g'); 
    hold on; 
end

subplot(4,2,z+1)
    plot(l:m:n, my_posterior(l:m:n), 'g'); 
    hold on;

%chain with IC3
load('chain_1_IC_3_final')
load('my_posterior_1_IC_3_final')

for i = 1:z
    subplot(4,2, i)
    plot(l:m:n, chain(l:m:n, i), 'b'); 
    hold on;
end

subplot(4,2,z+1)
    plot(l:m:n, my_posterior(l:m:n), 'b'); 
    hold on;


%chain with IC4
load('chain_1_IC_4_final')
load('my_posterior_1_IC_4_final')

for i = 1:z
    subplot(4,2, i)
    plot(l:m:n, chain(l:m:n, i), 'm'); 
end

subplot(4,2,z+1)
plot(l:m:n, my_posterior(l:m:n), 'm'); 
ylim([-2380, -2320])

   fig.Units = 'inches';
   fig.PaperPosition = [1 1 7 8.5];
   fig.Position = fig.PaperPosition;

   print(fig,'figS1.tiff','-dtiff','-r300')

end
