function Y = init_con_objective(initial_X0,t,model_params)
initial_conditions = [57.58e6-initial_X0(1)-initial_X0(2)-initial_X0(3)-initial_X0(4)-initial_X0(5);
    initial_X0;0];
[t_opt,x_opt] = ode45(@(t,initial_conditions) SEIQRDP_model(t,initial_conditions,model_params),t,initial_conditions);
tc_pred = x_opt(:,4)+x_opt(:,5)+x_opt(:,6);
Y = x_opt(:,4:6)';
end