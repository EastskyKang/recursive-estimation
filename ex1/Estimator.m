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
    
    posEst = [0,0];
    oriEst = 0;
    driftEst = 0;
    posVar = 0.333*estConst.TranslationStartBound^2*ones(1,2);
    oriVar = 0.333*estConst.RotationStartBound^2;
    driftVar = 0.333*estConst.GyroDriftStartBound^2;
    Pm = diag([posVar, oriVar, driftVar]);
    
    estState = updateEstState(estState, tm, posEst, oriEst, driftEst, posVar, oriVar, driftVar, Pm);

    return;
end


%% Mode 2: Estimator iteration.
% If we get this far tm is not equal to zero, and we are no longer
% initializing.  Run the estimator.
%Prior update
[xp_k, Pp_k] = priorUpdate(estState,actuate,sense,tm,estConst);

[xm, Pm] = measUpdate(estState,actuate,sense,tm,estConst, xp_k, Pp_k);

posEst = xm(1:2);
oriEst = xm(3);
driftEst = xm(4);
posVar = [Pm(1,1), Pm(2,2)];
oriVar = Pm(3,3);
driftVar = Pm(4,4);

estState = updateEstState(estState, tm, posEst, oriEst, driftEst, posVar, oriVar, driftVar, Pm);

end

function [xm, Pm] = measUpdate(estState,actuate,sense,tm,estConst, xp, Pp)

%Measurement update
z = [sense(2); sense(3); sense(1)]; %Measurement [zc, zg, zd]
H = [0, 0, 1, 0;
    0, 0, 1, 1;
    xp(1)/sqrt(xp(1)^2+xp(2)^2), xp(2)/sqrt(xp(1)^2+xp(2)^2), 0, 0];

M = eye(3);

R = diag([estConst.CompassNoise, estConst.GyroNoise, estConst.DistNoise]);

K = Pp * H'/(H * Pp * H'+M * R * M');
hk = h(xp, [0, 0, 0]);
err = K(:, z < inf) * (z(z < inf) - hk(z < inf));
xm = xp + err'; % update for only measured variable

Pm = (eye(4)-K(:, z < inf)*H(z < inf, :))*Pp;

end

function z = h(x, w)
% xp : x_hat p
% w : noise
z = [x(3) + w(1);
     x(3) + x(4) + w(2);
     sqrt(x(1)^2 + x(2)^2) + w(3)]; 
 
end

function [xk, Pk] = priorUpdate(estState,actuate,sense,tm,estConst)

tspan = [estState.tm, tm];

x0 = [estState.posEst, estState.oriEst, estState.driftEst]';
[tx,x] = ode45(@(t,x) q(x, estState, actuate, sense, tm, estConst), tspan, x0);

xk = x(end, :);

[~,P] = ode45(@(t,P) Pdot(t, tx, x, P, estState, actuate, sense , estConst), tspan, estState.Pm);

Pk = reshape(P(end, :), [4,4]);

end

function x_dot = q(x, estState,actuate, sense, tm, estConst)
% State update
sv = estConst.WheelRadius*actuate(1);
st = sv*cos(actuate(2));
sr = -sv*sin(actuate(2))/estConst.WheelBase;

x_dot = [st*cos(x(3)); st*sin(x(3)); sr; 0];

end

function P_dot = Pdot(t, tx, x, P, estState, actuate, sense , estConst)
% State update

x_hat = interp1(tx, x, t);
sv = estConst.WheelRadius*actuate(1);
st = sv*cos(actuate(2));

% Variance update
A = [0, 0, -st*sin(x_hat(3)), 0;
    0, 0, st*cos(x_hat(3)), 0;
    0, 0, 0, 0;
    0, 0, 0, 0];

L = [st*cos(x_hat(3)), 0;
    st*sin(x_hat(3)), 0;
    -sv*sin(actuate(2))/estConst.WheelBase, 0;
    0, 1];

Qc = estConst.VelocityInputPSD;
Qb = estConst.GyroDriftPSD;
Q = diag([Qc, Qb]);

P = reshape(P, 4, 4);
P_dot = A*P + P*A'+L*Q*L';
P_dot = P_dot(:);

end

function  [estState] = updateEstState(estState, tm, posEst, oriEst, driftEst, posVar, oriVar, driftVar, Pm)

    estState.tm = tm;
    estState.posEst = posEst;
    estState.oriEst = oriEst;
    estState.driftEst = driftEst;
    estState.posVar = posVar;
    estState.oriVar = oriVar;
    estState.driftVar = driftVar;
    estState.Pm = Pm;

end