function pos = multi_level_posterior_peak(parameters, x, k, indiv_PI_1, indiv_PI_2, indiv_PI_3, indiv_SI_DF_1, indiv_SI_DF_2, indiv_SI_DF_3, indiv_SI_DHF_1, indiv_SI_DHF_2, indiv_SI_DHF_3)

if k == 1 || k == 2 || k == 3
     x1 = cat(2, x(1:3),x(4),  0, 0, x(5), x(7), x(8), 10^-3.2, 1); %PI
     x2 = cat(2, x(1:3),x(4), parameters(13), x(6),x(5), x(7), x(8), 10^-3.2, 1); %SI DF
     x3 = cat(2, x(1:3),x(4), parameters(13), x(6), x(5), x(7), x(8), 10^-3.2, 1);  %SI DHF
elseif k == 4 
    beta_DHF = x(1) + x(9); 
    x1 = cat(2, x(1), x(2:4),  0, 0, x(5), x(7), x(8), 10^-3.5, 1); %PI
    x2 = cat(2, x(1), x(2:4), parameters(13), x(6),x(5), x(7), x(8), 10^-3.2, 1); %SI DF
    x3 = cat(2, beta_DHF, x(2:4), parameters(13), x(6), x(5), x(7), x(8), 10^-3.2, 1);  %SI DHF
elseif k == 5
     beta_1_DHF = x(1) + x(9); 
     beta_2_DHF = x(10); 
     beta_3_DHF = x(11); 
    x1 = cat(2, x(1), x(2:4),  0, 0, x(5), x(7), x(8), 10^-3.2, 1); %PI
    x2 = cat(2, x(1), x(2:4), parameters(13), x(6),x(5), x(7), x(8), 10^-3.2, 1); %SI DF
    x3 = cat(2, beta_1_DHF, x(2:4), parameters(13), x(6), x(5), beta_2_DHF + x(7), beta_3_DHF + x(8), 10^-3.2, 1);  %SI DHF
end
     
L_PI = multi_level_likelihood(parameters, x1, k, 100, indiv_PI_1, indiv_PI_2, indiv_PI_3, 9, 2, 2);
L_SI_DF = multi_level_likelihood(parameters, x2, k, 100, indiv_SI_DF_1, indiv_SI_DF_2, indiv_SI_DF_3, 29, 3, 2);
L_SI_DHF = multi_level_likelihood(parameters, x3, k, 100, indiv_SI_DHF_1, indiv_SI_DHF_2, indiv_SI_DHF_3, 10, 6, 2);

pos = L_PI + L_SI_DF + L_SI_DHF + priors(x);

end
