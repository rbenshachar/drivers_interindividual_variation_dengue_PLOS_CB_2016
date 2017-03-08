function void = plot_model_1_correlations(void)

%SUPPLEMENTARY FIGURE 2 in PLOS CB 2016 paper
%plots correlations of parameters of model 1

close all; 

load('chain_1_IC_1_final');
l = 150000; m = 100; n = 300000; 
my_chain = chain; 

fig = figure();

names = {'\beta', '\kappa', 'q','\sigma_I', 'q_T','log(V_0)'};
z = 6; 


x = 6; y = 6; 
temp = 1; 

for j = 1:z
    for i = 1:z        
        
       if i < j
            subplot(x,y,temp)
            plot(my_chain(l:m:n, i), my_chain(l:m:n, j), 'k.'); 
            set(gca, 'FontSize', 14);
            
           if i == 1
               ylabel(names(j))
           end
           
           if j == 6
               xlabel(names(i)) 
          end

       end
           temp = temp + 1; 
           

   end    
end

  fig.Units = 'inches';
   fig.PaperPosition = [1 1 13.5 8];
   fig.Position = fig.PaperPosition;

   print(fig,'figS2.tiff','-dtiff','-r300')
