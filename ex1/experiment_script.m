
num_iterations = 50;

simConst = SimulationConstants();
estConst = EstimatorConstants();
doplot=false;
seed = 0;

for i=1:num_iterations
    trackerror(i) = run(simConst,estConst,doplot,seed);
    fprintf('iteration: %d \n', i);
end

histogram(trackerror, 0:0.1:3);