%% This function takes three inputs
% x - a set of parameters
% t - the number of time-steps you wish to simulate

function f = slirdoutput(x,t,data)

% Here is a suggested framework for x.  However, you are free to deviate
% from this if you wish.

% set up transmission constants
susceptible_infection_rate = x(1);
mortality_rate = x(2);
recovery_rate = x(3);
vaccination_rate = x(4);
lockdown_infection_rate = x(5);

% set up initial conditions
ic_susceptible = x(6);
ic_lockdown = x(7);
ic_infected = x(8);
ic_recovered = x(9);
ic_dead = x(10);

% Set up SLIRD within-population transmission matrix
A = [1-susceptible_infection_rate-vaccination_rate, 0, 0, 0, 0;
    0, 1-lockdown_infection_rate, 0, 0, 0;
    susceptible_infection_rate, lockdown_infection_rate, 1-mortality_rate-recovery_rate, 0, 0;
    vaccination_rate, 0, recovery_rate, 1, 0;
    0, 0, mortality_rate, 0, 1];

% The next line creates a zero vector that will be used a few steps.
B = zeros(5,1);

% Set up the vector of initial conditions
x0 = [ic_susceptible; ic_lockdown; ic_infected; ic_recovered; ic_dead];

% Here is a compact way to simulate a linear dynamical system.
% Type 'help ss' and 'help lsim' to learn about how these functions work!!
sys_slird_base = ss(A,B,eye(5),zeros(5,1),1);
y = lsim(sys_slird_base,zeros(t,1),linspace(0,t-1,t),x0);
model = [y(:,3), y(:,5)];

% return a "cost".  This is the quantitity that you want your model to
% minimize.  Basically, this should encapsulate the difference between your
% modeled data and the true data. Norms and distances will be useful here.
% Hint: This is a central part of this case study!  choices here will have
% a big impact!
f = norm(model - data);

end