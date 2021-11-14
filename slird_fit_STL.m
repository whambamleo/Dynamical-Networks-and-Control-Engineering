
% Here is an example that reads in infection and fatalities from STL City
% and loads them into a new matrix covidstlcity_full
% In addition to this, you have other matrices for the other two regions in question
COVID_STL = COVID_MO(COVID_MO.name == "St. Louis",:);

% Change depending on what city we are modeling
start_range = 291;
end_range = 594;
COVIDdata = double(table2array(COVID_STL(:,(3:4))))./2805473;
COVIDdata = COVIDdata(start_range:end_range, :);
timesteps = table2array(COVID_STL(:,1));
timesteps = timesteps(start_range:end_range);
t = height(COVIDdata);

% The following line creates an 'anonymous' function that will return the cost (i.e., the model fitting error) given a set
% of parameters.  There are some technical reasons for setting this up in this way.
% Feel free to peruse the MATLAB help at
% https://www.mathworks.com/help/optim/ug/fmincon.html
% and see the sectiono on 'passing extra arguments'
% Basically, 'sirafun' is being set as the function siroutput (which you
% will be designing) but with t and coviddata specified.
sirafun= @(x)slirdoutput(x,t,COVIDdata);

%% set up rate and initial condition constraints
% Set A and b to impose a parameter inequality constraint of the form A*x < b
% Note that this is imposed element-wise
% If you don't want such a constraint, keep these matrices empty.
A = [];
b = [];

%% set up some fixed constraints
% Set Af and bf to impose a parameter constraint of the form Af*x = bf
% Hint: For example, the sum of the initial conditions should be
% constrained
% If you don't want such a constraint, keep these matrices empty.
Af = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
bf = 1;

%% set up upper and lower bound constraints
% Set upper and lower bounds on the parameters
% lb < x < ub
% here, the inequality is imposed element-wise
% If you don't want such a constraint, keep these matrices empty.
ub = [0.1, 0.1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1]';
lb = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0, 0, 0, 0, 0]';

% Specify some initial parameters for the optimizer to start from
% [susceptible_infection_rate, mortality_rate, recovery_rate, vaccination_rate, lockdown_infection_rate, initial_SLIRD]
x0 = [0.05, 0.01, 0.1, 0, 0, 0.9, 0, 0.1, 0, 0];

% This is the key line that tries to opimize your model parameters in order to
% fit the data
x = fmincon(sirafun,x0,A,b,Af,bf,lb,ub);
disp(x);

Y_fit = slirdmodel(x,t);

close all;

% Make some plots that illustrate your findings.
model_i = Y_fit(:,3);
model_d = Y_fit(:,5);
measured_i = COVIDdata(:,1);
measured_d = COVIDdata(:,2);

% Data
figure;
plot(timesteps, model_i, timesteps, model_d, timesteps, measured_i, timesteps, measured_d);
title('Fitted COVID Model vs. COVID Data - St. Louis', 'FontSize', 18);     % Title
xlabel('Date', 'FontSize', 18);              % Axis Labels
ylabel('Proportion of Population', 'FontSize', 18);
legend('Model I', 'Model D', 'Measured I', 'Measured D');

figure;
plot(Y_fit);
title('Fitted COVID Model - St. Louis', 'FontSize', 18);     % Title
xlabel('Days (since 2020-03-07)', 'FontSize', 18);              % Axis Labels
ylabel('Proportion of Population', 'FontSize', 18);
legend('Model S', 'Model L', 'Model I', 'Model R', 'Model D');