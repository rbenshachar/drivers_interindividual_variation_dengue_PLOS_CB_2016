%plot_serotype_viral_load_data

%FIGURE 1 OF PLOS CB PAPER

close all; 

x = 3; y = 3; 
%PI DF

[pat_num_1_PI_DF, time_1_PI_DF, data_1_PI_DF, ~] = viral_load_data(1,1,1);

[pat_num_2_PI_DF, time_2_PI_DF, data_2_PI_DF, ~] = viral_load_data(2,1,1);

[pat_num_3_PI_DF, time_3_PI_DF, data_3_PI_DF, ~] = viral_load_data(3,1,1);

fig = figure();

subplot(x, y, 1); 
for i = 1:length(pat_num_1_PI_DF)
           plot(time_1_PI_DF(i,:), data_1_PI_DF(i,:),'k','LineStyle', '-');
            hold on; 
end
plot(0:10, log10(1500).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
hold on; 
plot(0:10, log10(15000).*ones(11,1),  'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
xlim([0, 10]); ylim([2, 12]); 
ylabel({'DENV-1 viral load', '(log copies/ml)'}); 
title('1^° DF'); 
text(0.8,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)
text(0.65,0.7,'N = 15','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)

subplot(x, y, 4); 
for i = 1:length(pat_num_2_PI_DF)
           plot(time_2_PI_DF(i,:), data_2_PI_DF(i,:),'k','LineStyle', '-');
            hold on; 
end
plot(0:10, log10(1500).*ones(11,1),  'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
hold on; 
plot(0:10, log10(15000).*ones(11,1),  'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
xlim([0, 10]); ylim([2, 12]);   
ylabel({'DENV-2 viral load', '(log copies/ml)'});
text(0.8,0.9,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)
text(0.65,0.7,'N = 5','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)

subplot(x, y, 7); 
for i = 1:length(pat_num_3_PI_DF)
           plot(time_3_PI_DF(i,:), data_3_PI_DF(i,:),'k','LineStyle', '-');
            hold on; 
end
plot(0:10, log10(1500).*ones(11,1),  'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
hold on; 
plot(0:10, log10(15000).*ones(11,1),  'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
xlim([0, 10]); ylim([2, 12]); 
ylabel({'DENV-3 viral load', '(log copies/ml)'}); 
xlabel('days since onset of symptoms')
text(0.8,0.9,'(g)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)
text(0.65,0.7,'N = 6','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)

%SI DF

[pat_num_1_SI_DF, time_1_SI_DF, data_1_SI_DF, ~] = viral_load_data(1,2,1);

[pat_num_2_SI_DF, time_2_SI_DF, data_2_SI_DF, ~] = viral_load_data(2,2,1);

[pat_num_3_SI_DF, time_3_SI_DF, data_3_SI_DF, ~] = viral_load_data(3,2,1);

[pat_num_4_SI_DF, time_4_SI_DF, data_4_SI_DF, ~] = viral_load_data(4,2,1);

subplot(x, y, 2); 
for i = 1:length(pat_num_1_SI_DF)
           plot(time_1_SI_DF(i,:) , data_1_SI_DF(i,:),'k','LineStyle', '-');
            hold on; 
end
plot(0:10, log10(1500).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
hold on; 
plot(0:10, log10(15000).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
xlim([0, 10]); ylim([2, 12]); 
title('2^° DF')
text(0.8,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)
text(0.65,0.7,'N = 91','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)

subplot(x, y, 5); 
for i = 1:length(pat_num_2_SI_DF)
           plot(time_2_SI_DF(i,:), data_2_SI_DF(i,:),'k','LineStyle', '-');
            hold on; 
end
plot(0:10, log10(1500).*ones(11,1),'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
hold on; 
plot(0:10, log10(15000).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
xlim([0, 10]); ylim([2, 12]);     
text(0.8,0.9,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)
text(0.65,0.7,'N = 24','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)

subplot(x, y, 8); 
for i = 1:length(pat_num_3_SI_DF)
           plot(time_3_SI_DF(i,:) , data_3_SI_DF(i,:),'k','LineStyle', '-');
            hold on; 
end
plot(0:10, log10(1500).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
hold on; 
plot(0:10, log10(15000).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
xlim([0, 10]); ylim([2, 12]); 
xlabel('days since onset of symptoms'); 
text(0.8,0.9,'(h)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)
text(0.65,0.7,'N = 23','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)

%SI DHF
[pat_num_1_SI_DHF, time_1_SI_DHF, data_1_SI_DHF, ~] = viral_load_data(1,2,2);

[pat_num_2_SI_DHF, time_2_SI_DHF, data_2_SI_DHF, ~] = viral_load_data(2,2,2);

[pat_num_3_SI_DHF, time_3_SI_DHF, data_3_SI_DHF, ~] = viral_load_data(3,2,2);

subplot(x, y, 3); 
for i = 1:length(pat_num_1_SI_DHF)
           plot(time_1_SI_DHF(i,:), data_1_SI_DHF(i,:),'k','LineStyle', '-');
            hold on; 
end
plot(0:10, log10(1500).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
hold on; 
plot(0:10, log10(15000).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
xlim([0, 10]); ylim([2, 12]); 
text(0.8,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)
text(0.65,0.7,'N = 33','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)
title('2^° DHF')

subplot(x, y, 6); 
for i = 1:length(pat_num_2_SI_DHF)
           plot(time_2_SI_DHF(i,:) , data_2_SI_DHF(i,:),'k','LineStyle', '-');
            hold on; 
end
plot(0:10, log10(1500).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
hold on; 
plot(0:10, log10(15000).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
xlim([0, 10]); ylim([2, 12]);     
text(0.8,0.9,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)
text(0.65,0.7,'N = 21','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)


subplot(x, y, 9); 
for i = 1:length(pat_num_3_SI_DHF)
           plot(time_3_SI_DHF(i,:), data_3_SI_DHF(i,:),'k','LineStyle', '-');
            hold on; 
end
plot(0:10, log10(1500).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
hold on; 
plot(0:10, log10(15000).*ones(11,1), 'Color',  [0.5 0.5 0.5], 'LineStyle', '-.'); 
xlim([0, 10]); ylim([2, 12]); 
xlabel('days since onset of symptoms')
text(0.8,0.9,'(i)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)
text(0.65,0.7,'N = 10','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 12)

%FOR FORMATTING AND SAVING AS TIFF
fig.Units = 'inches';
fig.PaperPosition = [1 1 7.5 5];
fig.Position = fig.PaperPosition;

print(fig,'fig1.tiff','-dtiff','-r300')
