function [t,x] = seiqrdp_simulate_multiple_strain_waves(date_step, date_start, date_end, model_params,X0)
% Simulating the single strain SEIQRDP Covid model
X0 =    [X0(1), X0(1);
        1*X0(2),0*X0(2);
        0*X0(2),1*X0(2);
        0*X0(2),0*X0(2);
        1*X0(3),0*X0(3);
        0*X0(3),1*X0(3);
        0*X0(3),0*X0(3);
        1*X0(4),0*X0(4);
        0*X0(4),1*X0(4);
        0*X0(4),0*X0(4);
          X0(5),  X0(5);
          X0(6),  X0(6);
          X0(7),  X0(7)];
num_waves = size(date_start,1); % Get the number of waves by the number of start dates
date_range = datetime(date_start(1,:)):date_step:datetime(date_end(end,:));
tspan = 0:1:days(date_range(end) - date_range(1));
t_offset = [1;days(datetime(date_start(2,:)) - datetime(date_start(1,:)))];
x_par = zeros(length(tspan),13,num_waves);
t_par = zeros(length(tspan),num_waves);
for i = 1:num_waves
    tspan_wave = 0:1:length(tspan)-t_offset(i);
    wave_X0 = X0(:,i);
    [t_par(t_offset(i):end,i),x_par(t_offset(i):end,:,i)] = ode45(@(tspan_wave,wave_X0) SEIQRDP_multi_model(tspan_wave,wave_X0,model_params,i), ...
                                                        tspan_wave,wave_X0);
end
% x1_padded = [x_par(:,:,1);zeros(length(tspan)-t_offset(2),13)];
% x2_padded = [zeros(t_offset(2),13);x_par(:,:,2)];
x = x_par(:,:,1) + x_par(:,:,2);
% x = [x_par(:,:,1);x_par(:,:,2)];
t = datetime(date_start(1,:)):1:datetime(date_end(end,:));
% y = TODO add output vector