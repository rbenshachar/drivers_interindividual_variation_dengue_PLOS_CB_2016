function pos = posterior_CM_peak(parameters,x, k, indiv_PI, indiv_SI_DF, indiv_SI_DHF)

    if k == 1 %OAS     
        x1 = cat(2, x(1:3), x(4) , 0, 0, x(5), 10^-3.2, 1); %PI
        delta_DF = x(7) + parameters(13); 
        delta_DHF = parameters(13);  
        x2 = cat(2, x(1:3), x(4), delta_DF, x(6), x(5), 10^-3.2, 1); %SI DF
        x3 = cat(2, x(1:3), x(4), delta_DHF, x(6), x(5), 10^-3.2, 1); %SI DHF
    elseif k == 2 %ADE
        beta_DF = x(1); 
        beta_DHF = x(1) + x(7); 
        x1 = cat(2, x(1:3), x(4),  0, 0, x(5), 10^-3.2, 1); %PI
        x2 = cat(2, beta_DF, x(2:3), x(4), parameters(13), x(6), x(5), 10^-3.2, 1); %SI DF
        x3 = cat(2, beta_DHF, x(2:3), x(4), parameters(13), x(6), x(5),10^-3.2, 1); %SI DHF
    elseif k == 3
        beta_DHF = x(1) + x(7); 
        kappa_DHF = x(2); 
        kappa_DF = x(2) + x(8); 
        x1 = cat(2, x(1), kappa_DF, x(3), x(4),  0, 0, x(5), 10^-3.2, 1); %PI
        x2 = cat(2, x(1), kappa_DF, x(3), x(4), parameters(13), x(6), x(5), 10^-3.2, 1); %SI DF
        x3 = cat(2, beta_DHF, kappa_DHF, x(3), x(4), parameters(13), x(6), x(5),10^-3.2, 1); %SI DHF
    elseif k ==4
        beta_SI = x(1) + x(7); 
        q_p = x(3) + x(8); 
        x1 = cat(2, x(1:2), q_p, x(4),  0, 0, x(5), 1e-3); %PI
        x2 = cat(2, beta_SI, x(2), x(8), x(4), parameters(13), x(6), x(5), 1e-3); %SI
    end

    if k == 4
        indiv = indiv_PI;
         L_PI = likelihood_integrated(parameters, x1, 100, indiv, 9);
         L_SI = likelihood_integrated(parameters, x2, 100, indiv_SI_DHF, 28);
        pos = L_PI + L_SI + priors(x);
    else
       indiv = indiv_PI; 
        L_PI = likelihood_integrated(parameters, x1, 100, indiv, 13);
        indiv = indiv_SI_DF;  
        L_SI_DF = likelihood_integrated(parameters, x2, 100, indiv, 34);
        indiv = indiv_SI_DHF;  
        L_SI_DHF = likelihood_integrated(parameters, x3, 100, indiv, 18);
        pos = L_PI + L_SI_DF  + L_SI_DHF + priors(x);
    end
    
   

end