function [t,x] = seiqrdp_multi_simulate(date_step, date_start, date_end, model_params,X0)
% Simulating the single strain SEIQRDP Covid model
dt = date_step;
date_range = datetime(date_start):dt:datetime(date_end);
tspan = 0:1:days(date_range(end)-date_range(1));
X0 = [X0(1);
        1*X0(2);
        0*X0(2);
        0*X0(2);
        1*X0(3);
        0.05*X0(3);
        0.01*X0(3);
        1*X0(4);
        0.05*X0(4);
        0.01*X0(4);
        X0(5);
        X0(6);
        X0(7)];
[t,x] = ode45(@(tspan,X0) SEIQRDP_multi_model(tspan,X0,model_params,1),tspan,X0);
% y = TODO add output vector