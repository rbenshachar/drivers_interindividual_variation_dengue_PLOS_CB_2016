function pos = posterior(parameters,x, k, indiv_PI, indiv_SI)

%estimating beta, kappa, q, V0, deltaT, qT

    if k == 0       
        x1 = cat(2, x(1:3), 5.9, 0, 0, x(4), 10^(x(5)), 1); %PI
        x2 = x1; %SI
    elseif k == 1
        x1 = cat(2, x(1:3), 5.9, 0, 0,x(4), 10^(x(6)), 1); %PI
        x2 = cat(2, x(1:3), 5.9, parameters(13), x(5), x(4), 10^(x(6)), 1);
    elseif k == 2
        x1 = cat(2, x(1:3), 5.9, 0, 0,x(4), 10^(x(6)), 0.5); %PI
        x2 = cat(2, x(1:3), 5.9, parameters(13), x(5), x(4), 10^(x(6)), 0.5);
    elseif k == 3
        x1 = cat(2, x(1:3), 5.9, 0, 0,x(4), 10^(x(6)),2); %PI
        x2 = cat(2, x(1:3), 5.9, parameters(13), x(5), x(4), 10^(x(6)), 2); 
    elseif k == 4
        beta_SI = x(1) + x(6); 
        q_PI = x(2) + x(7); 
         x1 = cat(2, x(1),10, q_PI, x(3), 0, 0, x(4), 1e-3, 1); %PI
         x2 = cat(2, beta_SI, 10, x(2:3), parameters(13), x(5), x(4), 1e-3, 1);
    elseif k == 5
        beta_SI = x(1) + x(5); 
         x1 = cat(2, x(1),10, .11, x(2), 0, 0, x(4), 1e-3, 1); %PI
         x2 = cat(2, beta_SI, 10, .05, x(2), parameters(13), x(4), x(3), 1e-3, 1);
    end
    
    if k == 4 || k == 5
        L_PI = likelihood_integrated(parameters, x1, 100, indiv_PI, 11);
        L_SI = likelihood_integrated(parameters, x2, 100, indiv_SI, 90);
    else
        L_PI = likelihood_integrated(parameters, x1, 100, indiv_PI, 26);
        L_SI = likelihood_integrated(parameters, x2, 100, indiv_SI, 202);
    end

 
    if k == 4
         pos = L_PI + L_SI  + log(lognpdf(x(3), log(5.9), 0.15)); 
    elseif k == 5
        pos = L_PI + L_SI  + log(lognpdf(x(2), log(5.9), 0.15));
    else
        pos = L_PI + L_SI;
    end

end
