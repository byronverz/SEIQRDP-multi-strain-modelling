function [t,x] = seiqrdp_simulate(date_step, date_start, date_end, model_params,X0)
% Simulating the single strain SEIQRDP Covid model
dt = date_step;
date_range = datetime(date_start):dt:datetime(date_end);
tspan = 0:1:days(date_range(end)-date_range(1));

[t,x] = ode45(@(tspan,X0) SEIQRDP_model(tspan,X0,model_params,1),tspan,X0);
% y = TODO add output vector