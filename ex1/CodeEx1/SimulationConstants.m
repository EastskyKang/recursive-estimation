function const = SimulationConstants()
% const = UnknownConstants()
% 
% Define the constants used in the simulation.  These constants are not 
% accessible to the estimator.
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
% Raffaello D'Andrea, Michael , Lukas Hewing
% michaemu@ethz.ch
% lhewing@ethz.ch
%
% --
% Revision history
% [15.04.11, ST]    first version
% [07.05.13, MH]    first version
% [24.04.15, MM]    2015 version
% [05.05.17, LH]    2017 version


%% Speed

% The approximate maximum forward speed, in m/s.
const.MaxSpeedTranslation = 1.0;

% The approximate maximum rotational speed, in rad/s.
const.MaxSpeedRotation = 1.0;

% Minimum time for a segment, in seconds.
% Should be multiple of sampling time.
const.minSegTime = 0.5;

% Maximum time for a segment, in seconds.
% Should be multiple of sampling time.
const.maxSegTime = 2;


%% Robot kinematic constants

% The right and left wheel radius (W), in meter.
const.WheelRadius = 0.1;

% The wheel base (B), in meters.
const.WheelBase = 0.5;

%% Noise properties

% The compass sensor noise (w_c), normally distributed with zero mean 
% and variance \sigma_c^2, units rad^2.
const.CompassNoise = 0.1; % const.CompassNoise = \sigma_c^2

% The gyro sensor noise (w_g), normally distributed with zero mean 
% and variance \sigma_g^2, units rad^2.
const.GyroNoise = 0.005; % const.GyroNoise = \sigma_g^2

% The distance sensor noise (w_d), normally distributed with zero mean
% and variance \sigma_d^2, units rad^2.
const.DistNoise = 0.05; % const.DistNoise = \sigma_d^2

% Power spectral density of noise on wheel angular velocity commands (Q_v); 
% multiplicative noise, unit (rad/s)^2/Hz
const.VelocityInputPSD = 0.1; % const.VelocityInputPSD = Q_v

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

%% Times

% Number of samples of the simulation, the total duration of the simulation 
% is then const.N*const.sampleContinuous in seconds.
const.N = 500;

% The sample time for the continuous dynamics, in seconds.
const.sampleContinuous = 0.1;

% The min sample time for the compass, in seconds.
const.sampleCompassMin = 0.5;

% The max sample time for the compass, in seconds.
const.sampleCompassMax = 1.0;

% The min sample time for the gyro, in seconds.
const.sampleGyroMin = 0.1;

% The max sample time for the gyro, in seconds.
const.sampleGyroMax = 0.2;

% The min sample time for the position sensors, in seconds.
const.samplePosMin = 0.5;

% The max sample time for the position sensors, in seconds.
const.samplePosMax = 1.0;

