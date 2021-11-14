%% Clean Enviroment
clear all;
close all;

%% Changable Variables
% We assume that the following 4 parameters describe the dynamics of the
% situation
infection_rate = 0.05;  % Proportion of susceptible who become infected
mortality_rate = 0.01;  % Proportion of infected who die
true_recovery = 0.1;    % Proportion of infected who become recovered
part_recovery = 0.04;   % Proportion of infected who become susceptible

% We assume that in the begining, everyone is either infected or
% susceptible; no one has recovered nor died initially.
initial_infected = 0.1; % Proportion of infected at t = 0

% Time steps
num_steps = 1000;

%% Dynamics Calculations
% Dynamics Matrix
A = [1-infection_rate, part_recovery, 0, 0;
    infection_rate, 1-mortality_rate-true_recovery-part_recovery, 0, 0;
    0, true_recovery, 1, 0;
    0, mortality_rate, 0, 1];

% Time step vector
t = linspace(0, num_steps, num_steps + 1);

% Initialize output matrix with initial state
initial = [1-initial_infected; initial_infected; 0; 0];
output = [initial, zeros(4, num_steps)];

% Calculate each step
% Start at t = 1 since we already have initial state, by logical indexing
% this is the second index. End at t = timestep
for timestep=t(2:end)
    % The next timestep is the result of the previous timestep times the
    % dynamics matrix
    output(:,timestep + 1) = A * output(:,timestep);
end

%% Plotting
figure;                                             % Create a new figure
plot(t, output(1,:), t, output(2,:),...             % Plot all propportions
    t, output(3,:), t, output(4,:));
title('SIRD Model Simulation', 'FontSize', 18);     % Title
xlabel('Time (steps)', 'FontSize', 18);             % Axis Labels
ylabel('Proportion of Population', 'FontSize', 18);
legend('S', 'I', 'R', 'D', 'FontSize', 12);         % Legend