function [posEst,oriEst,driftEst, posVar,oriVar,driftVar,estState] = Estimator(estState,actuate,sense,tm,estConst)
% [posEst,oriEst,driftEst, posVar,oriVar,driftVar,estState] = 
%   Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function will be called in two different modes:
% If tm==0, the estimator is initialized; otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k), [1x2]-vector
%                   actuate(1): u_v, drive wheel angular velocity
%                   actuate(2): u_r, drive wheel angle
%   sense           sensor measurements z(k), [1x3]-vector, INF if no
%                   measurement
%                   sense(1): z_d, distance measurement
%                   sense(2): z_c, compass measurement
%                   sense(3): z_g, gyro measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConstants.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): x position estimate
%                   posEst(2): y position estimate
%   oriEst          orientation estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
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
% [19.04.11, ST]    first version by Sebastian Trimpe
% [30.04.12, PR]    adapted version for spring 2012, added unknown wheel
%                   radius
% [06.05.13, MH]    2013 version
% [23.04.15, MM]    2015 version
% [14.04.16, MM]    2016 version
% [05.05.17, LH]    2017 version


%% Mode 1: Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    pbar = estConst.TranslationStartBound;
    rbar = estConst.RotationStartBound;
    bbar = estConst.GyroDriftStartBound;      
    
    % mean
    % xhat_m[0] = x0
    xhat_m0 = [0, 0, 0, 0];
    
    % variance
    % Pm[0] = P0
    % variance for uniform zero mean RV = 1/3 x_bar ^2 
    vars = 1/3 * [pbar, pbar, rbar, bbar].^2;
    Pm0 = diag(vars);
    
    % independent
    diag_Pm0 = diag(Pm0);
    
    % Replace the following:
    posEst      = [xhat_m0(1), xhat_m0(2)];
    oriEst      = xhat_m0(3);
    driftEst    = xhat_m0(4);
    posVar      = [diag_Pm0(1), diag_Pm0(2)];
    oriVar      = diag_Pm0(3);
    driftVar    = diag_Pm0(4);
    
    % estimate state
    estState.x  = xhat_m0(1);
    estState.y  = xhat_m0(2);
    estState.r  = xhat_m0(3);
    estState.b  = xhat_m0(4);
    estState.Pm = Pm0;
    estState.tm = tm;

    return;
end


%% Mode 2: Estimator iteration.
% If we get this far tm is not equal to zero, and we are no longer
% initializing.  Run the estimator.

% -------------------------------------------------------------------------
% S1: prior update
% -------------------------------------------------------------------------

% mean
% dxhat(t) = q(xhat(t), 0, t)     
% tspan: (k-1)T <= t <= kT 
% initial cond: xhat((k-1)T) = xhat_m(k-1)

% tspan
tm_prev = estState.tm;          % (k-1)T
tspan = [tm_prev, tm];          % TODO CHECK

% initial cond
% xhat_m[k-1]
xhat_m_prev = [estState.x; estState.y; estState.r; estState.b];   

% ode
[t_xhat, xhat] = ode45(...
    @(t,xhat) q(t, xhat, actuate, [0; 0], estConst), ...
    tspan, xhat_m_prev);

% xhat_p[k] = xhat(kT)
xhat_p = xhat(end, :)';

% -------------------------------------------------------------------------
% variance
% dP(t) = A(t) P(t) + P(t) A'(t) + L(t) Q_c L'(t)
% tspan: (k-1)T <= t <= kT 
% initial cond: P((k-1)T) = Pm(k-1)

% tspan
tm_prev = estState.tm;
tspan = [tm_prev, tm];          % TODO CHECK

% initial cond
% Pm[k-1]
Pm_prev = estState.Pm;

% ode
[~, p] = ode45(...
    @(t, p) PmatODE(t, p, actuate, xhat, t_xhat, estConst), ...
    tspan, Pm_prev(:));

% Pp[k] = P(kT)
Pp = reshape(p(end, :)', 4, 4);

% -------------------------------------------------------------------------
% S2: measurement update
% -------------------------------------------------------------------------

% H matrix
denom = sqrt(xhat_p(1)^2 + xhat_p(2)^2);    % sqrt(x^2 + y^2)

H = [...
    0                   0                   1       0;
    0                   0                   1       1;
    xhat_p(1)/ denom    xhat_p(2) / denom   0       0];

% M matrix
M = eye(3);

% R matrix 
var_c = estConst.CompassNoise;
var_g = estConst.GyroNoise;
var_d = estConst.DistNoise;

R = diag([var_c, var_g, var_d]);

% -------------------------------------------------------------------------
% Kalman gain

joseph_form = true;

% z[k] 
z = [sense(2); sense(3); sense(1)];     % column vector (zc, zg, zd)'

% remove invalid measurements 
mask = ~isinf(z);
z = z(mask);

% H matrix invalid measurement removed
H = H(mask, :);

% M matrix invalid measurement removed
M = M(mask, :);

% K[k] = Pp[K] H'[k] (H[k] Pp[k] H'[k] + M[k] R M'[k])^-1
K = Pp * H' / (H * Pp * H' + M * R * M');

% -------------------------------------------------------------------------
% mean
% xhat_m[k] = xhat_p[k] + K[k] (z[k] - h_k(xhat_p[k], 0)) 

% h_k(xhat_p[k], 0) invalid measurement removed
h_xhat_p = h(xhat_p, [0; 0; 0]);
h_xhat_p = h_xhat_p(mask);

xhat_m = xhat_p + K * (z - h_xhat_p);

% variance

if joseph_form
    % Pm[k] = (I - K[k] H[k]) Pp[k] (I - K[k] H[k])' + K[k] M[k] R[k] M'[k] K'[k]
    Pm = (eye(4) - K * H) * Pp * (eye(4) - K * H)' + K * M * R * M' * K';
else
    % Pm[k] = (I - K[k] H[k]) Pp[k]
    Pm = (eye(4) - K * H) * Pp;
end

diag_Pm = diag(Pm);

% -------------------------------------------------------------------------
% Replace the following:
posEst      = [xhat_m(1), xhat_m(2)];
oriEst      = xhat_m(3);
driftEst    = xhat_m(4);
posVar      = [diag_Pm(1), diag_Pm(2)];
oriVar      = diag_Pm(3);
driftVar    = diag_Pm(4);

% estimate state
estState.x  = xhat_m(1);
estState.y  = xhat_m(2);
estState.r  = xhat_m(3);
estState.b  = xhat_m(4);
estState.Pm = Pm;
estState.tm = tm;

end

function [dstatedt] = q(t, state, actuate, proNoise, estConst)
    W = estConst.WheelRadius;
    B = estConst.WheelBase;
    
    uv = actuate(1);
    ur = actuate(2);
    
    vv = proNoise(1);
    vb = proNoise(2);
    
    x = state(1);
    y = state(2);
    r = state(3);
    b = state(4);
    
    dstatedt = zeros(4, 1);
    
    dstatedt(1) = W * uv * (1 + vv) * cos(ur) * cos(r);
    dstatedt(2) = W * uv * (1 + vv) * cos(ur) * sin(r);
    dstatedt(3) = - W / B * uv * (1 + vv) * sin(ur);
    dstatedt(4) = vb;
end

function [z] = h(state, meaNoise)
    x = state(1);
    y = state(2);
    r = state(3);
    b = state(4);
    
    wc = meaNoise(1);
    wg = meaNoise(2);
    wd = meaNoise(3);
    
    z = zeros(3, 1);
    
    z(1) = r + wc;
    z(2) = r + b + wg;
    z(3) = sqrt(x^2 + y^2) + wd;
end

function [dpdt] = PmatODE(t, p, actuate, xhat_v, t_xhat, estConst)
    % p is vectorized P 
    % xhat_v and t_xhat generate xhat(t) 
    
    W = estConst.WheelRadius;
    B = estConst.WheelBase;
    
    % xhat(t)
    xhat = interp1(t_xhat, xhat_v, t);
    rhat = xhat(3);
    
    % actuation
    uv = actuate(1);
    ur = actuate(2);
    
    % A matrix
    A = [...
        0   0   -W * uv * cos(ur) * sin( rhat )     0;
        0   0    W * uv * cos(ur) * cos( rhat )     0;
        0   0   0                                   0;
        0   0   0                                   0];
    
    % L matrix
    L = [...
        W * uv * cos(ur) * cos( rhat )      0;
        W * uv * cos(ur) * sin( rhat )      0;
        - W / B * uv * sin(ur)              0;
        0                                   1];
    
    % Q matrix
    Qv = estConst.VelocityInputPSD;
    Qb = estConst.GyroDriftPSD;
    Q = blkdiag(Qv, Qb);
        
    % reshape p to P
    P = reshape(p, 4, 4);

    % ode equation
    dPdt = A * P + P * A' + L * Q * L';
    
    % vectorize again
    dpdt = reshape(dPdt, size(p));
end
