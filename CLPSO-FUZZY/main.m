% Run CLPSOFUZZY once
projectDir = fileparts(mfilename('fullpath'));
addpath(genpath(projectDir));

objectiveFcn = @ackley;
iterations = 2000;
dims = 30;
bounds = [-32.768, 32.768];

[bestFit, bestInd, BestFitIter, BestIndIter, FAccessIter] = ...
    CLPSOFUZZY(objectiveFcn, iterations, dims, bounds);

fprintf('Best fitness: %.12g\n', bestFit);
disp('Best solution found:');
disp(bestInd);
fprintf('Function evaluations: %d\n', FAccessIter(end));

function fit = ackley(x)
    x = x(:).';
    d = numel(x);
    fit = -20 * exp(-0.2 * sqrt(sum(x.^2) / d)) ...
          - exp(sum(cos(2 * pi * x)) / d) + exp(1) + 20;
end
