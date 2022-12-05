function Population = int_pop(GenomeLength, ~, options)
% INT_POP Function that creates an initial population satisfying bounds and
% integer constraints

totalPopulation = sum(options.PopulationSize);

%IntCon constraints
IntCon = [1, 2,3,4];

range = options.PopInitRange;
lower = range(1,:);
span =  range(2,:) - lower;


Population = repmat(lower,totalPopulation,1 )+  ...
    repmat(span,totalPopulation,1) .* rand(totalPopulation, GenomeLength);

x = rand;
if x>=0.5
    Population(:,IntCon) = floor(Population(:, IntCon));
else
    Population(:,IntCon) = ceil(Population(:, IntCon));
end
Population = checkboundsIntGA(Population, range);