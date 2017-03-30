function pos = posterior_peak(parameters,x, k, indiv_PI, indiv_SI)

%estimating beta, kappa, deltaT, qT

    if k == 0       
        x1 = cat(2, x(1:3), x(4), 0, 0, x(5), 10^-3.2, 1); %PI
        x2 = x1; %SI
    elseif k == 1
         x1 = cat(2, x(1:3), x(4), 0, 0, x(5), 10^-3.2, 1); %PI
         x2 = cat(2, x(1:3), x(4),  parameters(13), x(6), x(5), 10^-3.2, 1);
    end

    L_PI = likelihood_integrated(parameters, x1, 100, indiv_PI, 13);
    L_SI = likelihood_integrated(parameters, x2, 100, indiv_SI, 52);

    pos = L_PI + L_SI + priors(x);

end
