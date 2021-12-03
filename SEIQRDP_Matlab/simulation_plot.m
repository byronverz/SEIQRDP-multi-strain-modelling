function simulation_plot(date_step, start_date,end_date,model_params)
% Plot simulation results of seiqrdp model
[x_opt,t] = seiqrdp_main(date_step,start_date,end_date,model_params);
subplot(1,2,1);
plot(t,x_opt(:,4:6)',date_range,ydata,'.');
legend("Quarantined fitted","Recovered fitted","Death fitted", ...
       "Quarantined data", "Recovered data", "Death data",...
    Location="best")