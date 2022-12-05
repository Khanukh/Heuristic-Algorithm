
FitnessFunction = @(x) vectorized_multiobjective(x);
nvars = length(ssinks);

lb = zeros(1,nvars) +1 ;
ub = zeros(1,nvars) +4 ;


rng default % for reproducibility
 ConstraintFunction = @simple_constraint;
 A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
intcon = zeros(1,nvars);
for i = 1 : nvars
    intcon(i)=i;
end

populationSize = 100;
generations = 30;
MaxStallGenerations = 10;
global Gvalue_num;
  global G_value;
   
    Gvalue_num =  1;

options = gaoptimset('PopulationSize',populationSize,...
        'Generations', generations ,...
        'TolFun', 170);

%   [x,fval] = gamultiobj(FitnessFunction,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);
tic 
[x,fval] = ga(FitnessFunction,nvars,A,b,Aeq,beq,lb,ub,ConstraintFunction,intcon,options);
toc
[x,fval]
plot(G_value)