function lambda = lambda_estimate(confirmed_data,recovered_data,initial_guess)
lambda_rate = diff(recovered_data)./confirmed_data(2:end);
figure("Name","Recovery rate function");
plot(lambda_rate);
lb = [0;0;0];
ub = [1;1;1];
t = 0:1:length(confirmed_data)-2;
opts = optimset('TolX',1e-15,'TolFun',1e-15,'Display','iter');
[lambda,resnorm, residuals] = lsqcurvefit(@(lambdas,t) lambda_func(t,lambdas), ...
                                         initial_guess,t,lambda_rate,lb,ub, ...
                                         opts);
% lambda = [0.113,0.104,0.007];
lambda
fprintf("Fitting resnorm: %d",resnorm);
figure("Name","Lambda Residuals");
plot(residuals');
figure("Name","Lambda fitted")
plot(t, lambda_rate,"r-.", ...
     t, lambda_func(t,lambda),"g-")
legend("Lambda data","Lambda fitted");
