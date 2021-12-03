function kappa = kappa_estimate(confirmed_data,deceased_data,initial_guess)
kappa_rate = diff(deceased_data)./confirmed_data(2:end);
figure("Name","Death rate function");
plot(kappa_rate);
lb = [0;0];
ub = [1;1];
t = 0:1:length(confirmed_data)-2;
opts = optimset('TolX',1e-15,'TolFun',1e-15,'Display','iter');
[kappa,resnorm, residuals] = lsqcurvefit(@(kappas,t) kappa_func(t,kappas), ...
                                         initial_guess,t,kappa_rate,lb,ub, ...
                                         opts);
% kappa = [0.002, 0.0];
kappa
fprintf("Fitting resnorm: %d",resnorm);
figure("Name","Kappa Residuals");
plot(residuals');
figure("Name","Kappa fitted")
plot(t, kappa_rate,"r-.", ...
     t, kappa_func(t,kappa),"g-")
legend("Kappa data","Kappa fitted");
