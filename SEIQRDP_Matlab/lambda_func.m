function Lambda = lambda_func(t,lambdas)
lambda1 = lambdas(1);
lambda2 = lambdas(2);
lambda3 = lambdas(3);
Lambda = lambda1 - lambda2*exp(-lambda3*t);
end