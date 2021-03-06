function trackErrorNorm=run(simConst,estConst,doplot,seed)
% run()
%
% Main function for Extended Kalman Filter programming exercise.
% It is used to simulate the true model, call the estimator and show the 
% results.
%
% Inputs:
%   simConst    constants used for the simulation (as in SimulationConstants.m)
%   estConst    estimator constants (as in EstimatorConstants.m), constants
%               used by the estimator
%   doplot      if doplot=true      plots are generated
%               if doplot=false     no plots are generated
%   seed        seed for the random number generator for reproducable results, 
%               if seed = 0 no seed is spedicfied
%
% Output:
%   trackErrorNorm  the tracking error e, as defined in the exercise sheet
%
% You can run the function simply with "run()" if you want to use the
% default parameters as provided by EstimatorConstants.m and
% SimulationConstants.m, and generate the plots.
%
%
%
% Class:
% Recursive Estimation
% Spring 2017
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Michael Muehlebach, Lukas Hewing
% michaemu@ethz.ch
% lhewing@ethz.ch
%
% --
% Revision history
% [15.04.11, ST]    first version by Sebastian Trimpe
% [30.04.12, PR]    2012 version, added unknown wheel radius
% [06.05.13, MH]    2013 version
% [24.04.15, MM]    2015 version
% [14.04.15, MM]    2016 version
% [05.05.17, LH]    2017 version

% clear command window, close figures
clc;
close all;

if nargin==0
   % Define the simulation constants that are used in simulation, but not 
   % accessible to the estimator.  
   simConst = SimulationConstants();
   
   % Define the constants accessible to the estimator.
   estConst = EstimatorConstants();
   
   % Generate plots by default.
   doplot=true;
   
   % use random seed
   seed = 0;
end



%% Setup

% Set the random number generator state.
% Uncomment to make results reproducable. This setting was used to generate
% the plot in the problem description.
if seed ~= 0
    rand('seed',seed);
    randn('seed',seed);
end


%% Simulation
% The function 'Simulator' simulates the robot kinematics and generates
% measurements.
[tm, loc, drift, input, sense] = Simulator( simConst );


%% Run the Estimator
N=simConst.N;

% Initialize the estimator.  
estState = [];
posEst = zeros(N,2);
oriEst = zeros(N,1);
driftEst = zeros(N,1);
posVar = zeros(N,2);
oriVar = zeros(N,1);
driftVar = zeros(N,1);
[posEst(1,:),oriEst(1),driftEst(1),posVar(1,:),oriVar(1),driftVar(1),estState] = ...
    Estimator(estState,zeros(1,2),zeros(1,3),0,estConst);

% Call the estimator for each time step.
for n = 2:N
    [posEst(n,:),oriEst(n),driftEst(n),posVar(n,:),oriVar(n),driftVar(n),estState] = ...
        Estimator(estState,input(n,:),sense(n,:),tm(n),estConst);
end



%% The results
% Plots of the results.

% Calculate the total tracking error.
% Replace the following:
trackErrorNorm = sqrt(mean(sum((posEst - loc(:, 1:2)).^2, 2)));

if doplot
    % Add your plots here to debug the estimator and verify your
    % implementation.
    
    % 2d map (x, y)
    figure('Name', 'trajectory')
    plot(loc(:, 1), loc(:, 2), 'r.')
    hold on 
    plot(posEst(:, 1), posEst(:, 2), 'b.')
    plot(loc(1, 1), loc(1, 2), 'kX')
    hold off
    title('Y-X')
    legend('real', 'est', 'start pt')
    
    % orientation (vs time)
    figure('Name', 'orientation')
    plot(loc(:, 3))
    hold on 
    plot(oriEst)
    hold off
    title('Orientation')
    legend('real', 'est')
    
    % gyro drift (vs time)
    figure('Name', 'gyro')
    plot(drift)
    hold on 
    plot(driftEst)
    hold off
    title('Gyro drift')
    legend('real', 'est')
    
    % variances
    figure('Name', 'variances')
    plot(posVar(:, 1))
    hold on 
    plot(posVar(:, 2))
    plot(oriVar)
    plot(driftVar)
    hold off 
    title('Variance')
    legend('x', 'y', 'ori', 'drift')
end
    
return;

    