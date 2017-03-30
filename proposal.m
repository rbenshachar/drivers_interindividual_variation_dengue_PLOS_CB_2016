function p = proposal(means, my_cov)
        current_vals = means; 

        p = mvnrnd(current_vals, my_cov);
end

