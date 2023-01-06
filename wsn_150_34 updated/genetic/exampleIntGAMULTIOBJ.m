%% GAMULTIOBJ with integer constraints
%
% Example on how to use GAMULTIOBJ for mixed-integer programming with bound
% constraints

%% Problem setup
fitnessFunction = @srn;  % Function handle to the fitness function
numberOfVariables = 2;   % Number of decision variables
populationSize = 50;
stallGenLimit = 200;
generations = 500;

% Bound Constraints
lb = [-20, -20];       % Lower bound
ub = [20, 20];         % Upper bound
Bound = [lb; ub];

%% Solve the problem without integer constraints
options = gaoptimset('PopulationSize',populationSize,...
    'PopInitRange',Bound, 'StallGenLimit', stallGenLimit,...
    'Generations', generations);

[x, f, exitflag, output, population, score] = gamultiobj(fitnessFunction,...
    numberOfVariables, [], [], [], [], lb, ub, options);

%% Solve the problem with integer constraints
%Add one of the following lines to each of the custom functions:
%IntCon = [];       - no integer constraints
%IntCon = 1;        - x(1) of the solution will be integer
%IntCon = 2;        - x(2) of the solution will be integer
%IntCon = [1, 2];   - both cordinates of the solution will be integer

options = gaoptimset('PopulationSize',populationSize,...
    'CreationFcn', @int_pop,...
    'MutationFcn', @int_mutation,...
    'CrossoverFcn',@int_crossoverarithmetic,...
    'StallGenLimit', stallGenLimit,...
    'Generations', generations,...
    'PopulationSize',populationSize,...
    'PopInitRange', Bound);

[x1, f1, exitflag1, output1, population1, score1] = gamultiobj(fitnessFunction,...
     numberOfVariables, [], [], [], [], lb, ub, options);
 
%% Plot the results

% Compare the Pareto fronts for the problems with and without integer
% constraints
fh1 = figure;
scatter(f(:, 1), f(:, 2), 'fill');
hold all
scatter(f1(:, 1), f1(:, 2), 'fill');
xlabel('f_1');
ylabel('f_2');
legend({'only bounds', 'bounds+integer'});
title('Pareto front comparison with and without integer constraints');

% Compare the solutions of the problems with and without integer
% constraints
fh2 = figure('units', 'normalized', 'position',[0.2, 0.2, 0.6, 0.6] );
scatter(x(:, 1), x(:, 2), 'fill');
hold all
scatter(x1(:, 1), x1(:, 2), 'fill');
xlabel('x');
ylabel('y');
xlim([-20, 20]);
ylim([-20, 20]);
set(gca, 'Xtick', -20:2:20, 'Ytick', -20:2:20, 'XminorTick', 'on', 'YminorTick', 'on');
grid on
legend({'only bounds', 'bounds+integer'});
title('Solution comparison with and without integer constraints');