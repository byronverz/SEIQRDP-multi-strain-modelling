% Main script for the SEIQRDP model fitting and simulation
%% Section 1: Data extraction
% Extracting the confirmed, recovered and death data from downloaded data
import_opts = detectImportOptions("..\CovidModelling\data\time_series_covid19_confirmed_global.csv", ...
                           "ReadVariableNames",true, ...
                           "NumHeaderLines",0, ...
                           "VariableNamingRule","preserve", ...
                           "TextType","string");
confirmed_table = readtable("..\CovidModelling\data\time_series_covid19_confirmed_global.csv",import_opts);
sa = confirmed_table(confirmed_table.("Country/Region") == "South Africa",:);
south_africa_confirmed = removevars(sa,1:4);
sac = south_africa_confirmed.Variables;
recovered_table = readtable("..\CovidModelling\data\time_series_covid19_recovered_global.csv",import_opts);
sa = recovered_table(recovered_table.("Country/Region") == "South Africa",:);
south_africa_recovered= removevars(sa,1:4);
sar = south_africa_recovered.Variables;
death_table = readtable("..\CovidModelling\data\time_series_covid19_deaths_global.csv",import_opts);
sa = death_table(death_table.("Country/Region") == "South Africa",:);
south_africa_dead = removevars(sa,1:4);
sad = south_africa_dead.Variables;
% Setting the date range for the data to use
start_date = datetime([2020,4,18]); % 18 April 2020
end_date = datetime([2020,10,1]); % 1 Oct 2020
offset_date = datetime([2020,1,22]); % Need this since the data starts on 22 Jan
south_africa_data_start = days(start_date- offset_date);
south_africa_data_end = days(end_date - offset_date);
total_confirmed = sac(south_africa_data_start:south_africa_data_end);
recovered = sar(south_africa_data_start:south_africa_data_end);
deaths = sad(south_africa_data_start:south_africa_data_end);
confirmed = total_confirmed-(recovered+deaths);
date_range = start_date:1:end_date;
figure();
subplot(2,1,1);
hold on
grid on
plot(date_range,total_confirmed,'Color',[0 0.4470 0.7410],'LineStyle','none','Marker','.');
plot(date_range,recovered,'Color',[0.4940 0.1840 0.5560],'LineStyle','none','Marker','.');
plot(date_range,deaths,'Color',[0 0 0],'LineStyle','none','Marker','.');
plot(date_range,confirmed,'Color',[0.8500 0.3250 0.0980],'LineStyle','none','Marker','.');
legend("Total","Recovered","Deaths","Quarantined","Location","northwest");
hold off
subplot(2,1,2);
hold on
grid on
plot(date_range(1:end-1),diff(total_confirmed),"b-", ...
     date_range(1:end-1),diff(recovered),"g-", ...
     date_range(1:end-1),diff(deaths),"k-");
legend("Daily new confirmed","Daily new recovered","Daily new deaths","Location","northwest");
hold off
%% Section 2: Parameter fitting
% Fitting the nonlinear, time-variant SEIQRDP model for single strain case
% Setting the initial conditions
N = 57.58e6; % Total population
Q0 = 1500; % Initial quarantined (can also be seen as initial confirmed infectious cases)
% Q0 = 0;
E0 = 700; % Initial exposed
% E0 = 0;
I0 = 2200; % Initial infectious
% I0 = 0;
R0 = 500; % Initial recovered
D0 = 90; % Initial deceased
P0 = 0; % Initial insusceptible
S0 = N - E0 - D0 - P0 - Q0; % Initial susceptible
X0 = [S0;
      E0;
      I0;
      Q0;
      R0;
      D0;
      P0];
% We need to estimate the time dependant parameters first to increase the
% estimation of the final parameter fitting

% Fitting the mortality rate (kappa) which can be seen as the rate of new
% deceased cases for every case moving out of quarantine
kappa = kappa_estimate(confirmed,deaths,[0.002,0.0]);

% Fitting the recovery rate (lambda) which can be seen as the rate of new
% recovered cases for every case moving out of quarantine
lambda = lambda_estimate(confirmed,recovered,[0.113,0.104,0.007]);
params0 = [0.01,...       % Beta1
           0.6,...       % Beta2
           0.01,...       % Beta3
           0.2, ...    % Gamma
           0.18];     % Delta
set_params = [10e-6,...     % Alpha
              lambda(1),... % Lambda1
              lambda(2),... % Lambda1
              lambda(3),... % Lambda1
              kappa(1),...  % Kappa1
              kappa(2)];    % Kappa2

%%
fitted_params = fit_seiqrdp_model(confirmed,recovered,deaths,X0,params0,set_params);
fprintf("Alpha %d \n",set_params(1))
fprintf("Beta1 %2.3f \n",fitted_params(1))
fprintf("Beta2 %2.3f \n",fitted_params(2))
fprintf("Beta3 %2.3f \n",fitted_params(3))
fprintf("Gamma %2.3f \n",fitted_params(4))
fprintf("Delta %2.3f \n",fitted_params(5))
fprintf("Lambda1 %2.3f \n",set_params(2))
fprintf("Lambda2 %2.3f \n",set_params(3))
fprintf("Lambda3 %2.3f \n",set_params(4))
fprintf("Kappa1 %2.3f \n",set_params(5))
fprintf("Kappa2 %2.3f \n",set_params(6))
fitted_parameters = [set_params(1),fitted_params,set_params(2:end)];
%% Section 3: Simulation

[t_fitted,x_fitted] = seiqrdp_simulate(1,start_date,end_date,fitted_parameters,X0);
total_confirmed_fitted = x_fitted(:,4)+x_fitted(:,5)+x_fitted(:,6);
figure("Name","Fitting results");
hold on
grid on
plot(date_range,total_confirmed_fitted,'Color',[0 0.4470 0.7410]); 
plot(date_range,x_fitted(:,4),'Color',[0.8500 0.3250 0.0980]);
plot(date_range,x_fitted(:,5),'Color',[0.4940 0.1840 0.5560]);
plot(date_range,x_fitted(:,6),'Color',[0 0 0]);
plot(date_range,total_confirmed,'Color',[0 0.4470 0.7410],'LineStyle','none','Marker','.')
plot(date_range,confirmed,'Color',[0.8500 0.3250 0.0980], 'LineStyle','none','Marker','.');
plot(date_range,recovered,'Color',[0.4940 0.1840 0.5560], 'LineStyle','none','Marker','.');
plot(date_range,deaths,'Color',[0 0 0],'LineStyle','none','Marker','.')
legend("Total confirmed","Confirmed","Recovered","Deaths","Location","best")
hold off
%%
ydata = [confirmed;recovered;deaths];
SStot = sum((ydata-mean(ydata)).^2);                            % Total Sum-Of-Squares
SSres = sum((ydata-x_fitted(:,4:6)').^2);                         % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot
%%
new_date_end = datetime([2021,3,1]);
new_data_end = days(new_date_end - offset_date);
new_date_range = start_date:1:new_date_end;
total_new_confirmed = sac(south_africa_data_start:new_data_end);
new_recovered = sar(south_africa_data_start:new_data_end);
new_death = sad(south_africa_data_start:new_data_end);
new_confirmed = total_new_confirmed - (new_recovered+new_death);
[~,x_fitted_new] = seiqrdp_simulate(1,start_date,new_date_end,fitted_parameters,X0);
total_confirmed_fitted_new = x_fitted_new(:,4)+x_fitted_new(:,5)+x_fitted_new(:,6);
figure("Name","Extended Forecast");
hold on
grid on
plot(new_date_range,total_confirmed_fitted_new,'Color',[0 0.4470 0.7410])
plot(new_date_range,x_fitted_new(:,4),'Color',[0.8500 0.3250 0.0980]);
plot(new_date_range,x_fitted_new(:,5),'Color',[0.4940 0.1840 0.5560]);
plot(new_date_range,x_fitted_new(:,6),'Color',[0 0 0]);
plot(new_date_range,total_new_confirmed,'Color',[0 0.4470 0.7410],'LineStyle','none','Marker','.')
plot(new_date_range,new_confirmed,'Color',[0.8500 0.3250 0.0980], 'LineStyle','none','Marker','.');
plot(new_date_range,new_recovered,'Color',[0.4940 0.1840 0.5560], 'LineStyle','none','Marker','.');
plot(new_date_range,new_death,'Color',[0 0 0],'LineStyle','none','Marker','.');
legend("Total","Confirmed","Recovered","Deaths",'Location','northwest')
hold off
t_range = 0:1:(south_africa_data_end-south_africa_data_start);
beta_range = zeros(1,length(t_range));
kappa_range = zeros(1,length(t_range));
lambda_range = zeros(1,length(t_range));
for ti = t_range
    beta_range(ti+1) = beta_func(ti,fitted_parameters(2),fitted_parameters(3),fitted_parameters(4),1);
    kappa_range(ti+1) = kappa_func(ti, fitted_parameters(10:end));
    lambda_range(ti+1) = lambda_func(ti,fitted_parameters(7:9));
end

fig = figure("Name","Time variant parameters");
set(fig,'defaultAxesColorOrder',[[0 0.0 0.7410];[0 0 0]]);
hold on
yyaxis left
plot(date_range,beta_range,'b-');
plot(date_range,lambda_range,'b-.');
yyaxis right
plot(date_range,kappa_range,'k-');
hold off
legend("Spreading rate ($\beta$)", ...
    "Mortality rate ($\kappa$)", ...
    "Recovery rate ($\lambda$)",'Interpreter','latex');
%% Section 4: Multiple Waves
% We will use the first wave model parameters and simulate for multiple
% waves of one strain
wave_start_dates = [2020,4,18;2020,10,8]; % Two start dates for 2 waves
wave_end_dates = [2020,10,8;2021,3,30]; % 2 end dates for each wave
[t_multi_wave, x_multi_wave,beta_multi] = seiqrdp_simulate_multiple_waves(1,wave_start_dates,wave_end_dates, ...
                                                                fitted_parameters,X0);
multi_wave_total_sim = x_multi_wave(:,4)+x_multi_wave(:,5)+x_multi_wave(:,6);
multi_wave_data_end = days(datetime(wave_end_dates(end,:)) - offset_date);
multi_date_range = start_date:1:new_date_end;
multi_wave_total_data = sac(south_africa_data_start:multi_wave_data_end);
recovered_multi_wave = sar(south_africa_data_start:multi_wave_data_end);
death_multi_wave = sad(south_africa_data_start:multi_wave_data_end);
confirmed_multi_wave = multi_wave_total_data - (recovered_multi_wave+death_multi_wave);
figure("Name","Multiwave Simulation")
hold on
plot(t_multi_wave,multi_wave_total_sim,'Color',[0 0.4470 0.7410]);
plot(t_multi_wave,x_multi_wave(:,4),'Color',[0.8500 0.3250 0.0980]);
plot(t_multi_wave,x_multi_wave(:,5),'Color',[0.4940 0.1840 0.5560]);
plot(t_multi_wave,x_multi_wave(:,6),'Color',[0 0 0]);
plot(t_multi_wave,multi_wave_total_data,'Color',[0 0.4470 0.7410],'LineStyle','none','Marker','.')
plot(t_multi_wave,confirmed_multi_wave,'Color',[0.8500 0.3250 0.0980], 'LineStyle','none','Marker','.');
plot(t_multi_wave,recovered_multi_wave,'Color',[0.4940 0.1840 0.5560], 'LineStyle','none','Marker','.');
plot(t_multi_wave,death_multi_wave,'Color',[0 0 0],'LineStyle','none','Marker','.');
legend("Total","Confirmed","Recovered","Deaths",'Location','northwest')
hold off
figure("Name","Time varying parameters")
plot(t_multi_wave,beta_multi)
legend("$\beta$",'Interpreter','latex')
%% Section 5: Multi-strain model data
% First fetch the normalized strain contributions
import_opts = detectImportOptions("../CovidModelling/normalized_variant_freq.csv");
normed_strain_freqs_table = readtable("../CovidModelling/normalized_variant_freq.csv",import_opts);
origin_col = '20A-D(Origin Clade)';
normed_strain_freqs_table.(origin_col) = [1-sum(normed_strain_freqs_table{1:23,2:end},2);zeros(16,1)];
variant_date_range = normed_strain_freqs_table.Var1;
normed_strain_freqs = table2array(normed_strain_freqs_table(:,2:end));
normed_strain_freqs = smoothdata(normed_strain_freqs);
figure("Name","Strain Frequency")
area(variant_date_range,normed_strain_freqs)
legend(normed_strain_freqs_table.Properties.VariableNames(2:end),'Location','northwest')
% Simplified frequency
import_opts = detectImportOptions("../CovidModelling/normalized_variant_freq_scaled.csv");
normed_strain_freqs_table_scaled = readtable("../CovidModelling/normalized_variant_freq_scaled.csv",import_opts);
origin_col = '20A-D(Origin Clade)';
normed_strain_freqs_table_scaled.(origin_col) = 1-sum(normed_strain_freqs_table_scaled{:,{'Beta','Delta'}},2);
variant_date_range = normed_strain_freqs_table_scaled.Var1;
normed_strain_freqs_scaled = table2array(normed_strain_freqs_table_scaled(1:24,{'Beta','20A-D(Origin Clade)','Delta'}));
normed_strain_freqs_scaled = smoothdata(normed_strain_freqs_scaled,1,'sgolay',5);
figure("Name","Simplified Strain Frequency")
area(variant_date_range(1:24),normed_strain_freqs_scaled)
colours = [0.8500 0.3250 0.0980;0.6350 0.0780 0.1840;0.3010 0.7450 0.9330];
colororder(colours)
legend('Beta','20A-D(Origin Clade)','Delta','Location','northwest')
%%
% Split out a fitting data wave into "variants" based on frequency data
end_ix = find((variant_date_range.Year==end_date.Year)&(variant_date_range.Month==end_date.Month),1);
variant_start_date = variant_date_range(1);
variant_end_date = variant_date_range(end_ix);
split_date_range = variant_start_date:1:variant_end_date;
variant_split_data = normed_strain_freqs_scaled(1:end_ix ,:);
variant_split_expanded = zeros(length(split_date_range),7);
for i = 1:1:3
    variant_split_expanded(:,i) = interp1(variant_date_range(1:end_ix),variant_split_data(:,i),split_date_range,'spline');
end
plot(split_date_range,variant_split_expanded,variant_date_range(1:end_ix),variant_split_data,'og')
%%
split_start_data = days(variant_start_date - offset_date);
split_end_data = split_start_data + length(split_date_range)-1;
split_total_confirmed_cases = sac(split_start_data:split_end_data);
split_data = zeros(length(split_date_range),6);
for i = 1:1:7
    split_data(:,i) = split_total_confirmed_cases'.*variant_split_expanded(:,i);
end
plot(split_date_range,split_data,split_date_range,split_total_confirmed_cases)
%%
s1_f = split_data(:,1)'./split_total_confirmed_cases;
s2_f = split_data(:,2)'./split_total_confirmed_cases;
area(split_date_range,[s2_f',s1_f'])
%% Section 5: Multi-Strain Multi-Wave model simulation
strain_start_dates = [2020,4,18;2020,9,15]; % Two start dates for 2 waves
strain_end_date = [2021,3,30]; % 1 end date for the end of the whole sim
[t_multi_strain_wave, x_multi_strain_wave] = seiqrdp_simulate_multiple_strain_waves(1,wave_start_dates,wave_end_dates, ...
                                                                fitted_parameters,X0);
confirmed_all_strains = x_multi_strain_wave(:,8) + x_multi_strain_wave(:,9) + x_multi_strain_wave(:,10);
figure("Name","Multi Model Simulation")
plot(t_multi_strain_wave ,x_multi_strain_wave(:,8:12)');
legend("Strain1:Confirmed","Strain2:Confirmed","Strain3:Confirmed", ...
       "Recovered","Deceased",'Location','best');
figure("Name","Multi Model Confirmed")
plot(t_multi_strain_wave,x_multi_strain_wave(:,8:10)')
legend("Strain1:Confirmed","Strain2:Confirmed","Strain3:Confirmed","Location","northwest")
multi_wave_strain_data_end = days(datetime(strain_end_date) - offset_date);
multi_date_range = start_date:1:new_date_end;
multi_wave_strain_total_data = sac(south_africa_data_start:multi_wave_strain_data_end);
recovered_multi_wave_strain = sar(south_africa_data_start:multi_wave_strain_data_end);
death_multi_wave_strain = sad(south_africa_data_start:multi_wave_strain_data_end);
confirmed_multi_wave_strain = multi_wave_strain_total_data - (recovered_multi_wave_strain+death_multi_wave_strain);
multi_wave_strain_total = confirmed_all_strains + x_multi_strain_wave(:,11) + x_multi_strain_wave(:,12);  
figure("Name","Multi-wave Multi-strain Simulation")
hold on
plot(t_multi_strain_wave,multi_wave_strain_total,'Color',[0 0.4470 0.7410]);
plot(t_multi_strain_wave,x_multi_strain_wave(:,8:10),'Color',[0.8500 0.3250 0.0980]);
plot(t_multi_strain_wave,x_multi_strain_wave(:,11),'Color',[0.4940 0.1840 0.5560]);
plot(t_multi_strain_wave,x_multi_strain_wave(:,12),'Color',[0 0 0]);
plot(t_multi_strain_wave,multi_wave_total_data,'Color',[0 0.4470 0.7410],'LineStyle','none','Marker','.')
plot(t_multi_strain_wave,confirmed_multi_wave,'Color',[0.8500 0.3250 0.0980], 'LineStyle','none','Marker','.');
plot(t_multi_strain_wave,recovered_multi_wave,'Color',[0.4940 0.1840 0.5560], 'LineStyle','none','Marker','.');
plot(t_multi_strain_wave,death_multi_wave,'Color',[0 0 0],'LineStyle','none','Marker','.');
legend("Total","Confirmed S_1","Confirmed S_2","Confirmed S_3", ...
    "Recovered","Deaths",'Location','best')
hold off
%%
no_var_states = x_multi_strain_wave(:,8:10)';
[t_multi_strain_wave, x_multi_strain_wave_var] = seiqrdp_simulate_multiple_strain_waves(1,wave_start_dates,wave_end_dates, ...
                                                                fitted_parameters,X0);
plus_ten_var_states = x_multi_strain_wave_var(:,8:10)';
%%
[t_multi_strain_wave, x_multi_strain_wave_var2] = seiqrdp_simulate_multiple_strain_waves(1,wave_start_dates,wave_end_dates, ...
                                                                fitted_parameters,X0);
minus_ten_var_states = x_multi_strain_wave_var2(:,8:10)';
%%
plot(t_multi_strain_wave,no_var_states,'-', ...
    t_multi_strain_wave,plus_ten_var_states,'--', ...
    t_multi_strain_wave,minus_ten_var_states,'--')
%%
% figure("Name","Time varying parameters")
% plot(t_multi_strain_wave,beta_multi)
% legend("$\beta$",'Interpreter','latex')
%%
% Multi strain frequency plot
multi_freq = [x_multi_strain_wave(:,9)./confirmed_all_strains,x_multi_strain_wave(:,8)./confirmed_all_strains,x_multi_strain_wave(:,10)./confirmed_all_strains];
figure("Name","Multi strain simulation frequency")
area(t_multi_strain_wave,multi_freq)
colours = [0.8500 0.3250 0.0980;0.6350 0.0780 0.1840];
colororder(colours)
legend("Strain2","Strain1")
%%
ydata = [confirmed_multi_wave;recovered_multi_wave;death_multi_wave];
xdata = [confirmed_all_strains,x_multi_strain_wave(:,11),x_multi_strain_wave(:,12)];
SStot = sum((ydata-mean(ydata)).^2);                            % Total Sum-Of-Squares
SSres = sum((ydata-xdata').^2);                         % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot





