function model_params = fit_seiqrdp_model(confirmed_data,recovered_data,death_data,x0,params0,set_params)
lower_bounds = [0.0,...     % Beta1
                0.0,...     % Beta2
                0.0,...     % Beta3
                0.2,...   % Gamma
                0.1];   % Delta
upper_bounds = [1,...      % Beta1
                1,...      % Beta2
                1,...      % Beta3
                1,...         % Gamma
                0.2];         % Delta

tc_ydata = confirmed_data+recovered_data+death_data;
% y_data = [tc_ydata;confirmed_data;recovered_data;death_data];
y_data = [confirmed_data;recovered_data;death_data];
t = 0:1:length(confirmed_data)-1;
optim_options = optimoptions('lsqcurvefit','Display','iter', ...
                             'MaxFunctionEvaluations',1e6, ...
                             'MaxIterations',1.0e3, ...
                             'FunctionTolerance',1e-6, ...
                             'StepTolerance',1e-6, ...
                             'OptimalityTolerance',1e-6);
[model_params,resnorm,residuals] = lsqcurvefit(@(params0,t) objective_func(params0,t,x0,set_params), ...
                                               params0,t,y_data,lower_bounds,upper_bounds, ...
                                               optim_options);
fprintf("Resnorm: %d \n",resnorm);
figure("Name","Model fitting residuals");
plot(residuals');
legend("Confirmed","Recovered","Deaths")
title("Parameter fitting residuals")