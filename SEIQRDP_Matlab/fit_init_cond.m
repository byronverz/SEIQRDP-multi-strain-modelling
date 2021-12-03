function X0 = fit_init_cond(model_params, initial_X0,confirmed,recoverd,deaths)
lower_bounds = [200,...               % E0
                100,...               % I0
                300,...              % Q0
                0,...              % R0
                0];                % D0
upper_bounds = [3000,...              % E0
                4000,...              % I0
                4000,...              % Q0
                1000,...              % R0
                1000];                % D0

tc_ydata = confirmed+recoverd+deaths;
% y_data = [tc_ydata;confirmed_data;recovered_data;death_data];
y_data = [confirmed;recoverd;deaths];
t = 0:1:length(confirmed)-1;
optim_options = optimoptions('lsqcurvefit','Display','off', ...
                             'MaxFunctionEvaluations',1e9, ...
                             'MaxIterations',1.25e3, ...
                             'FunctionTolerance',1e-9, ...
                             'StepTolerance',1e-9, ...
                             'OptimalityTolerance',1e-9,'Algorithm','trust-region-reflective');
[X0,resnorm,residuals] = lsqcurvefit(@(initial_X0,t) init_con_objective(initial_X0,t,model_params), ...
                                               initial_X0,t,y_data,lower_bounds,upper_bounds, ...
                                               optim_options);
% fprintf("Resnorm: %d \n",resnorm);
% figure("Name","Initial condition fitting residuals");
% plot(residuals');
% legend("Confirmed","Recovered","Deaths")
% title("Initial condition fitting residuals")