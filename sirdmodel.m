%% This function takes three inputs
% x - a set of parameters
% t - the number of time-steps you wish to simulate
% data - actual data that you are attempting to fit

function f = sirdmodel(x,t)

% set up transmission constants
infection_rate = x(1);
mortality_rate = x(2);
recovery_rate = x(3);

% set up initial conditions
ic_susceptible = x(4);
ic_infected = x(5);
ic_recovered = x(6);
ic_dead = x(7);

% Set up SIRD within-population transmission matrix
A = [1-infection_rate, 0, 0, 0;
    infection_rate, 1-mortality_rate-recovery_rate, 0, 0;
    0, recovery_rate, 1, 0;
    0, mortality_rate, 0, 1];

% The next line creates a zero vector that will be used a few steps.
B = zeros(4,1);

% Set up the vector of initial conditions
x0 = [ic_susceptible; ic_infected; ic_recovered; ic_dead];

% simulate the SIRD model for t time-steps
sys_sir_base = ss(A,B,eye(4),zeros(4,1),1);
y = lsim(sys_sir_base,zeros(t,1),linspace(0,t-1,t),x0);

f = y;

end