%% Paths

addpath('problems/cannon');
addpath('results_manager');
addpath('propagators');
addpath('stopping_criteria');
addpath('timestep_control');

%% Configuration

% Initial conditions
% vector estado inicial
w0 = [0; % Pos X
    0;   % Pos Y
    5;   % Vel X
    100  % Vel Y
    ];

% Constant acceleration
a_x = 0;
a_y = -9.81;

% Our problem function should be a function of just the state vector and
% the current time. cannon is a function of also the accelerations. So
% we are using here an anonymous function (or lambda) to create a function
% that always calls cannon passing w and t as variables and using constant
% values for a_x and a_y, the ones defined earlier.
problem = @(w, t) cannon(w, t, a_x, a_y);
% Lambda function. Recibe los argumentos w, t y devuelve el resultado de cannon(w, t, a_x, a_y)
% Se guarda en la constante problem

% Set the initial time and the current state vector to the initial one.
t0 = 0;
w = w0;

% The timestep calculator should be a function of just the state vector and
% the current time. Using constant_dt with a lambda function to get a
% constant timestep of 0.01 s.
dt_fun = @(w, t) constant_dt(w, t, 0.01);

% Try two time integrators: an explicit Euler method and an implicit Euler
% method.
expl_propagator = @fw_euler;
impl_propagator = @bw_euler;

% Set the stopping criteria.
condition = @stop_at_floor;

% Set the results manager. It has to be a function of the state vector, the
% current time and the previous results structure. We are using here a
% constant sampling period sampler, with a sampling period of 0.01 s.
manager = @(w, t, old_results) sample_results(w, t, old_results, 0.01);


%% Solve with an explicit method.
expl_results = solver(w0, t0, problem, expl_propagator, ...
    dt_fun, condition, manager);

%% Solve with an implicit method.
impl_results = solver(w0, t0, problem, impl_propagator, ...
    dt_fun, condition, manager);

%% Results postprocessing.
plot(expl_results(:, 2), expl_results(:, 3), 'r-');
hold on;
plot(impl_results(:, 2), impl_results(:, 3), 'b-');
xlabel('$x$ [m]', 'Interpreter', 'latex');
ylabel('$y$ [m]', 'Interpreter', 'latex');
legend('Explicit', 'Implicit');
hold off;