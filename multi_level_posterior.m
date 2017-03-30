function pos = multi_level_posterior(parameters, x, k, indiv_PI_1, indiv_PI_2, indiv_PI_3, indiv_SI_DF_1, indiv_SI_DF_2, indiv_SI_DF_3, indiv_SI_DHF_1, indiv_SI_DHF_2, indiv_SI_DHF_3)

if k == 4 %SS_beta + ADE
     beta_DHF = x(1) + x(8); 
     x1 = cat(2, x(1:3),5.9,  0, 0, x(4), x(6), x(7), 10^(x(9))); %PI
     x2 = cat(2, x(1:3),5.9, parameters(13), x(5),x(4), x(6), x(7),10^(x(9))); %SI DF
     x3 = cat(2, beta_DHF, x(2:3),5.9, parameters(13), x(5), x(4), x(6), x(7),10^(x(9)));  %SI DHF
elseif k == 5 %ADE for each serotype
     beta_1_DHF = x(1) + x(8); 
     beta_2_DHF = x(9); 
     beta_3_DHF = x(10); 
     
     x1 = cat(2, x(1:3),5.9,  0, 0, x(4), x(6), x(7), 10^(x(11))); %PI
     x2 = cat(2, x(1:3),5.9, parameters(13), x(5),x(4), x(6), x(7),10^(x(11))); %SI DF
     x3 = cat(2, beta_1_DHF, x(2:3),5.9, parameters(13), x(5), x(4),beta_2_DHF + x(6),beta_3_DHF + x(7),10^(x(11)));  %SI DHF
else %serotype-specific differences in beta, q, or qT
     x1 = cat(2, x(1:3),5.9,  0, 0, x(4), x(6), x(7), 10^(x(8))); %PI
     x2 = cat(2, x(1:3),5.9, parameters(13), x(5),x(4), x(6), x(7),10^(x(8))); %SI DF
     x3 = cat(2, x(1:3),5.9, parameters(13), x(5), x(4), x(6), x(7),10^(x(8)));  %SI DHF
end

    L_PI = multi_level_likelihood(parameters, x1, k, 100, indiv_PI_1, indiv_PI_2, indiv_PI_3, 15, 5, 6);
    L_SI_DF = multi_level_likelihood(parameters, x2, k, 100, indiv_SI_DF_1, indiv_SI_DF_2, indiv_SI_DF_3, 91, 24, 23);
    L_SI_DHF = multi_level_likelihood(parameters, x3, k, 100, indiv_SI_DHF_1, indiv_SI_DHF_2, indiv_SI_DHF_3, 33, 21, 10);

    pos = L_PI + L_SI_DF + L_SI_DHF; 

end
