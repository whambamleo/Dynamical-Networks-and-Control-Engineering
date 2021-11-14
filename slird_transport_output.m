%% This function takes three inputs
% x - a set of parameters
% t - the number of time-steps you wish to simulate

function f = slird_transport_output(x,t,data)

% Here is a suggested framework for x.  However, you are free to deviate
% from this if you wish.

% STL parameters
% set up transmission constants
stl_susceptible_infection_rate = x(1);
stl_mortality_rate = x(2);
stl_recovery_rate = x(3);
stl_vaccination_rate = x(4);
stl_lockdown_infection_rate = x(5);

% set up initial conditions
stl_ic_susceptible = x(6);
stl_ic_lockdown = x(7);
stl_ic_infected = x(8);
stl_ic_recovered = x(9);
stl_ic_dead = x(10);

% set up travel constants
stl_to_kc = x(11);
stl_to_mo = x(12);

% KC parameters
% set up transmission constants
kc_susceptible_infection_rate = x(13);
kc_mortality_rate = x(14);
kc_recovery_rate = x(15);
kc_vaccination_rate = x(16);
kc_lockdown_infection_rate = x(17);

% set up initial conditions
kc_ic_susceptible = x(18);
kc_ic_lockdown = x(19);
kc_ic_infected = x(20);
kc_ic_recovered = x(21);
kc_ic_dead = x(22);

% set up travel constants
kc_to_stl = x(23);
kc_to_mo = x(24);

% MO parameters
% set up transmission constants
mo_susceptible_infection_rate = x(25);
mo_mortality_rate = x(26);
mo_recovery_rate = x(27);
mo_vaccination_rate = x(28);
mo_lockdown_infection_rate = x(29);

% set up initial conditions
mo_ic_susceptible = x(30);
mo_ic_lockdown = x(31);
mo_ic_infected = x(32);
mo_ic_recovered = x(33);
mo_ic_dead = x(34);

% set up travel constants
mo_to_stl = x(35);
mo_to_mo = x(36);

% Set up SLIRD within-population transmission matrices
A_STL = [1-stl_susceptible_infection_rate-stl_vaccination_rate, 0, 0, 0, 0;
    0, 1-stl_lockdown_infection_rate, 0, 0, 0;
    stl_susceptible_infection_rate, stl_lockdown_infection_rate, 1-stl_mortality_rate-stl_recovery_rate, 0, 0;
    stl_vaccination_rate, 0, stl_recovery_rate, 1, 0;
    0, 0, stl_mortality_rate, 0, 1];
A_KC = [1-kc_susceptible_infection_rate-kc_vaccination_rate, 0, 0, 0, 0;
    0, 1-kc_lockdown_infection_rate, 0, 0, 0;
    kc_susceptible_infection_rate, kc_lockdown_infection_rate, 1-kc_mortality_rate-kc_recovery_rate, 0, 0;
    kc_vaccination_rate, 0, kc_recovery_rate, 1, 0;
    0, 0, kc_mortality_rate, 0, 1];
A_MO = [1-mo_susceptible_infection_rate-mo_vaccination_rate, 0, 0, 0, 0;
    0, 1-mo_lockdown_infection_rate, 0, 0, 0;
    mo_susceptible_infection_rate, mo_lockdown_infection_rate, 1-mo_mortality_rate-mo_recovery_rate, 0, 0;
    mo_vaccination_rate, 0, mo_recovery_rate, 1, 0;
    0, 0, mo_mortality_rate, 0, 1];
% Transport matrices
KC_to_STL = zeros(5, 5);
KC_to_STL(1, 1) = kc_to_stl;

MO_to_STL = zeros(5, 5);
MO_to_STL(1, 1) = mo_to_stl;

STL_to_KC = zeros(5, 5);
STL_to_KC(1, 1) = stl_to_kc;

MO_to_KC = zeros(5, 5);
MO_to_KC(1, 1) = mo_to_kc;

STL_to_MO = zeros(5, 5);
STL_to_MO(1, 1) = stl_to_mo;

KC_to_MO = zeros(5, 5);
KC_to_MO(1, 1) = kc_to_mo;

% Full matrix
A = [A_STL, KC_to_STL, MO_to_STL;
    STL_to_KC, A_KC, MO_to_KC;
    STL_to_MO, KC_to_MO, A_MO];

% The next line creates a zero vector that will be used a few steps.
B = zeros(15,1);

% Set up the vector of initial conditions
x0 = [stl_ic_susceptible; stl_ic_lockdown; stl_ic_infected; stl_ic_recovered; stl_ic_dead;...
    kc_ic_susceptible; kc_ic_lockdown; kc_ic_infected; kc_ic_recovered; kc_ic_dead;...
    mo_ic_susceptible; mo_ic_lockdown; mo_ic_infected; mo_ic_recovered; mo_ic_dead];

% Here is a compact way to simulate a linear dynamical system.
% Type 'help ss' and 'help lsim' to learn about how these functions work!!
sys_slird_base = ss(A,B,eye(5),zeros(5,1),1);
y = lsim(sys_slird_base,zeros(t,1),linspace(0,t-1,t),x0);
model = [y(:,3), y(:,5), y(:,8), y(:,10), y(:,13), y(:,15)];

% return a "cost".  This is the quantitity that you want your model to
% minimize.  Basically, this should encapsulate the difference between your
% modeled data and the true data. Norms and distances will be useful here.
% Hint: This is a central part of this case study!  choices here will have
% a big impact!
f = norm(model - data);

end