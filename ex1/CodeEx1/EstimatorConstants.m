function const = EstimatorConstants()
% const = KnownConstants()
% 
% Define the physical constants that are available to the estimator.
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
% [30.04.12, PR]    adapted version 2012, added unknown wheel radius
% [06.05.13, MH]    adapted version 2013
% [24.04.15, MM]    version 2015
% [14.04.16, MM]    version 2016
% [05.05.17, LH]    version 2017


%% Robot kinematic constants

% The right and left wheel radius (W_0), in meter.
const.WheelRadius = 0.1;

% The wheel base (B), in meters.
const.WheelBase = 0.5;

%% Noise properties

% The compass sensor noise (w_r), normally distributed with zero mean 
% and variance \sigma_r^2, units rad^2.
const.CompassNoise = 0.1; % const.CompassNoise = \sigma_r^2

% The gyro sensor noise (w_g), normally distributed with zero mean 
% and variance \sigma_g^2, units rad^2.
const.GyroNoise = 0.005; % const.GyroNoise = \sigma_g^2

% The distance sensor noise (w_d), normally distributed with zero mean
% and variance \sigma_d^2, units rad^2.
const.DistNoise = 0.05; % const.DistNoise = \sigma_d^2

% Power spectral density of noise on wheel angular velocity commands (Q_v); 
% multiplicative noise, unit (rad/s)^2/Hz
const.VelocityInputPSD = 0.1; % const.VelocityInputNoise = Q_v

% Power spectral density of Gyro drift (Q_b), units (rad/s)^2/Hz.
const.GyroDriftPSD = 0.1; % const.GyroDriftPSD = Q_b

%% Starting point

% The robot nominally starts at the origin, uniformly distributed with the
% following bound (\bar{p}), in meters.
const.TranslationStartBound = 0.5; % const.TranslationStartBound = \bar{p}

% The nominal orientation is also 0, and has a uniform distribution with
% the following bound (\bar{r}), in rad.
const.RotationStartBound = pi/8; % const.RotationStartBound = \bar{r}

% The initial gyro drift is exactly 0
const.GyroDriftStartBound = 0;