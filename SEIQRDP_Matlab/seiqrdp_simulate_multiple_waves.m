function [t,x,beta] = seiqrdp_simulate_multiple_waves(date_step, date_start, date_end, model_params,X0)
% Simulating the single strain SEIQRDP Covid model
dt = date_step;
num_waves = size(date_start,1); % Get the number of waves by the number of start dates
date_range = datetime(date_start(1,:)):dt:datetime(date_end(end,:));
tspan = 0:1:days(datetime(date_end(1,:)) - datetime(date_start(1,:)));
% t_offset = [1;days(datetime(date_start(2,:)) - datetime(date_start(1,:)))];
x_par = zeros(length(tspan),7,num_waves);
t_par = zeros(length(tspan),num_waves);
beta_range = zeros(num_waves,length(tspan));
for i = 1:num_waves
    tspan_wave = 0:1:days(datetime(date_end(i,:))-datetime(date_start(i,:)));
    [t_par(:,i),x_par(:,:,i)] = ode45(@(tspan_wave,X0) SEIQRDP_model(tspan_wave,X0,model_params,i), ...
                                                        tspan_wave,X0);
    X0 = x_par(end,:,i);
    for ti = tspan_wave
        beta_range(i,ti+1) = beta_func(ti,model_params(2),model_params(3),model_params(4),i);
    end
end
% x1_padded = [x_par(:,:,1);zeros(t_offset,7)];
% x2_padded = [zeros(t_offset,7);x_par(:,:,2)];
% x = x_par(:,:,1) + x_par(:,:,2);
x = [x_par(:,:,1);x_par(2:end,:,2)];
t = date_range;
beta= [beta_range(1,:),beta_range(2,2:end)];
% y = TODO add output vector