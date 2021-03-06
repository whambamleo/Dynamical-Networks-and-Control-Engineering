
% Here is an example that reads in infection and fatalities from STL City
% and loads them into a new matrix covidstlcity_full
% In addition to this, you have other matrices for the other two regions in question
COVID_STL = COVID_MO(COVID_MO.name == "St. Louis",:);
COVID_KC = ;
COVID_MO = ;

Total_Pop = ;

% Change depending on what city we are modeling
start_range = 291;
end_range = 594;
COVIDdataSTL = double(table2array(COVID_STL(:,(3:4))))./Total_Pop;
COVIDdataSTL = COVIDdataSTL(start_range:end_range, :);
COVIDdataKC = double(table2array(COVID_KC(:,(3:4))))./Total_Pop;
COVIDdataKC = COVIDdataKC(start_range:end_range, :);
COVIDdataMO = double(table2array(COVID_MO(:,(3:4))))./Total_Pop;
COVIDdataMO = COVIDdataMO(start_range:end_range, :);
COVIDdata = [COVIDdataSTL, COVIDdataKC, COVIDdataMO];
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
sirafun= @(x)slird_transport_output(x,t,COVIDdata);

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
Af = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0,...
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0,...
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0];
bf = 1;

%% set up upper and lower bound constraints
% Set upper and lower bounds on the parameters
% lb < x < ub
% here, the inequality is imposed element-wise
% If you don't want such a constraint, keep these matrices empty.
ub = [0.1, 0.1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1,...
    0.1, 0.1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1,...
    0.1, 0.1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1]';
lb = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0, 0, 0, 0, 0, 0, 0,...
    0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0, 0, 0, 0, 0, 0, 0,...
    0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0, 0, 0, 0, 0, 0, 0]';

% Specify some initial parameters for the optimizer to start from
% [susceptible_infection_rate, mortality_rate, recovery_rate, vaccination_rate, lockdown_infection_rate, initial_SLIRD]
x0 = [0.05, 0.01, 0.1, 0, 0, 0.9, 0, 0.1, 0, 0, 0, 0,...
    0.05, 0.01, 0.1, 0, 0, 0.9, 0, 0.1, 0, 0, 0, 0,...
    0.05, 0.01, 0.1, 0, 0, 0.9, 0, 0.1, 0, 0, 0, 0];

% This is the key line that tries to opimize your model parameters in order to
% fit the data
x = fmincon(sirafun,x0,A,b,Af,bf,lb,ub);
disp(x);

Y_fit = slird_transport_model(x,t);

close all;

% Make some plots that illustrate your findings.
model_i_STL = Y_fit(:,3);
model_d_STL = Y_fit(:,5);
model_i_KC = Y_fit(:,8);
model_d_KC = Y_fit(:,10);
model_i_MO = Y_fit(:,8);
model_d_MO = Y_fit(:,10);
measured_i_STL = COVIDdata(:,1);
measured_d_STL = COVIDdata(:,2);
measured_i_KC = COVIDdata(:,3);
measured_d_KC = COVIDdata(:,4);
measured_i_MO = COVIDdata(:,5);
measured_d_MO = COVIDdata(:,6);

% Data
figure;
plot(timesteps, model_i_STL, timesteps, model_d_STL,...
    timesteps, model_i_KC, timesteps, model_d_KC,...
    timesteps, model_i_MO, timesteps, model_d_MO,...
    timesteps, measured_i_STL, timesteps, measured_d_STL,...
    timesteps, measured_i_KC, timesteps, measured_d_KC,...
    timesteps, measured_i_MO, timesteps, measured_d_MO);
title('Fitted COVID Model vs. COVID Data - Transport', 'FontSize', 18);     % Title
xlabel('Date', 'FontSize', 18);              % Axis Labels
ylabel('Proportion of Population', 'FontSize', 18);
legend('STL Model I', 'STL Model D', 'KC Model I', 'KC Model D',...
    'MO Model I', 'MO Model D', 'STL Measured I', 'STL Measured D',...
    'KC Measured I', 'KC Measured D', 'MO Measured I', 'MO Measured D');

figure;
plot(Y_fit);
title('Fitted COVID Model - Transport', 'FontSize', 18);     % Title
xlabel('Days (since 2020-03-07)', 'FontSize', 18);              % Axis Labels
ylabel('Proportion of Population', 'FontSize', 18);
legend('STL Model S', 'STL Model L', 'STL Model I', 'STL Model R', 'STL Model D',...
    'KC Model S', 'KC Model L', 'KC Model I', 'KC Model R', 'KC Model D',...
    'MO Model S', 'MO Model L', 'MO Model I', 'MO Model R', 'MO Model D');