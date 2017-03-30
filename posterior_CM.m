function pos = posterior_CM(parameters,x, k, indiv_PI, indiv_SI_DF, indiv_SI_DHF)

    if k == 1 %OAS     
        x1 = cat(2, x(1:3), 5.9, 0, 0, x(4), 10^(x(7)), 1); %PI
        delta_DF = x(6) + parameters(13); 
        delta_DHF = parameters(13);  
        x2 = cat(2, x(1:3), 5.9, delta_DF, x(5), x(4), 10^(x(7)), 1); %SI DF
        x3 = cat(2, x(1:3), 5.9, delta_DHF, x(5), x(4), 10^(x(7)), 1); %SI DHF
    elseif k == 2 %ADE
        beta_DHF = x(1) + x(6); 
        x1 = cat(2, x(1:3), 5.9, 0, 0,x(4), 10^(x(7)), 1); %PI
        x2 = cat(2, x(1:3), 5.9, parameters(13), x(5), x(4), 10^(x(7)), 1); %SI DF
        x3 = cat(2, beta_DHF, x(2:3), 5.9, parameters(13), x(5), x(4), 10^(x(7)), 1); %SI DHF
    elseif k == 3
        x1 = cat(2, x(1:3), 5.9, 0, 0, x(4), 10^(x(8)), 1); %PI
        delta_DF = x(6) + parameters(13); 
        delta_DHF = parameters(13);  
        qT_DF = x(5); 
        qT_DHF = x(5) + x(7); 
        x2 = cat(2, x(1:3), 5.9, delta_DF, qT_DF, x(4), 10^(x(8)), 1); %SI DF
        x3 = cat(2, x(1:3), 5.9, delta_DHF, qT_DHF, x(4), 10^(x(8)), 1); %SI DHF
    end
    
    indiv = indiv_PI;       
    L_PI = likelihood_integrated(parameters, x1, 100, indiv, 26);
    indiv = indiv_SI_DF;  
    L_SI_DF = likelihood_integrated(parameters, x2, 100, indiv, 138);
    indiv = indiv_SI_DHF;  
    L_SI_DHF = likelihood_integrated(parameters, x3, 100, indiv, 64);
    
    pos = L_PI + L_SI_DF  + L_SI_DHF;

end