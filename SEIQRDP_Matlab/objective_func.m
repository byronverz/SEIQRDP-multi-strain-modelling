function Y = objective_func(params, t,X0,set_params)
model_parameters = [set_params(1),params,set_params(2:end)];

[t_opt,x_opt] = ode45(@(t,X0) SEIQRDP_model(t,X0,model_parameters,1),t,X0);
tc_pred = x_opt(:,4)+x_opt(:,5)+x_opt(:,6);
Y = x_opt(:,4:6)';
end