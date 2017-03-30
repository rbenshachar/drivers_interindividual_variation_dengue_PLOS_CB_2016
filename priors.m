function p = priors(proposal)

p = log(lognpdf(proposal(1,4), log(5.9), 0.15)); 

end
