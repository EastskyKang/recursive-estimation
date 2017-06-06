function [postParticles] = Estimator(prevPostParticles, sens, act, init)
% [postParticles] = Estimator(prevPostParticles, sens, act, init)
%
% The estimator function. The function will be called in two different
% modes: If init==1, the estimator is initialized. If init == 0, the
% estimator does an iteration for a single sample time interval Ts (KC.ts)
% using the previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurements and control inputs.
%
% You must edit this function.
%
% Inputs:
%   prevPostParticles   previous posterior particles at discrete time k-1,
%                       which corresponds to continuous time t = (k-1)*Ts
%                       The variable is a struct whose fields are arrays
%                       that correspond to the posterior particle states.
%                       The fields are: (N is number of particles)
%                       .x = (2xN) array with the x-locations (metres)
%                       .y = (2xN) array with the y-locations (metres)
%                       .h = (2xN) array with the headings (radians)
%                       The first row in the arrays corresponds to robot A.
%                       The second row corresponds to robot B.
%
%   sens                Sensor measurements at discrete time k (t = k*Ts),
%                       [4x1]-array, an Inf entry indicates no measurement
%                       of the corresponding sensor.
%                       sens(1): distance reported by sensor 1 (metres)hel
%                       sens(2): distance reported by sensor 2 (metres)
%                       sens(3): distance reported by sensor 3 (metres)
%                       sens(4): distance reported by sensor 4 (metres)
%
%   act                 Control inputs u at discrete time k-1, which are
%                       constant during a time interval Ts:
%                       u(t) = u(k-1) for (k-1)*Ts <= t < k*Ts
%                       [2x1]-array:
%                       act(1): velocity of robot A, u_A(k-1) (metres/second)
%                       act(2): velocity of robot B, u_B(k-1) (metres/second)
%
%   init                Boolean variable indicating wheter the estimator
%                       should be initialized (init = 1) or if a regular
%                       estimator update should be performed (init = 0).
%                       OPTIONAL ARGUMENT. By default, init = 0.
%
% Outputs:
%   postParticles       Posterior particles at discrete time k, which
%                       corresponds to the continuous time t = k*Ts.
%                       The variable is a struct whose fields are arrays
%                       that correspond to the posterior particle states.
%                       The fields are: (N is number of particles)
%                       .x = (2xN) array with the x-locations (metres)
%                       .y = (2xN) array with the y-locations (metres)
%                       .h = (2xN) array with the headings (radians)
%                       The first row in the arrays corresponds to robot A.
%                       The second row corresponds to robot B.
%
% Class:
% Recursive Estimation
% Spring 2017
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Michael Muehlebach
% michaemu@ethz.ch

% Check if init argument was passed to estimator:
if(nargin < 4)
    % if not, set to default value:
    init = 0;
end

%% Mode 1: Initialization
% Set number of particles:
N = 10; % obviously, you will need more particles than 10.
if (init)
    % Do the initialization of your estimator here!
    % These particles are the posterior particles at discrete time k = 0
    % which will be fed into your estimator again at k = 1
        % Do the initialization of your estimator here!
    % These particles are the posterior particles at discrete time k = 0
    % which will be fed into your estimator again at k = 1
    
    mask = rand(2,N) > 0.5;
    x0 = [2 * KC.L *ones(1, N); zeros(1, N)]; % Robots can be on ether side of the room
    y0 = KC.L* (mask); % Robots can be on ether side of the room
    h0 = [rand(1,N)*pi()/2+pi()/2;
        rand(1,N)*pi()/2];
    
    h0(1, mask(1, :)) = h0(1, mask(1, :)) + pi()/2;
    
    h0(2, mask(2, :)) = h0(2, mask(2, :)) + 3*pi()/2;
        
        
    postParticles.x = x0;
    postParticles.y = y0;
    postParticles.h = h0;
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If init = 0, we perform a regular update of the estimator.

% Implement your estimator here!

[xp, yp, hp] = priorUpdate(prevPostParticles, act);

[xm, ym, hm] = measUpdate(xp, yp, hp, sens);
%Roughening
[xm, ym, hm] = roughening(xm, ym, hm);

% Replace the following:
postParticles.x = xm;
postParticles.y = ym;
postParticles.h = hm;

end % end estimator

function [xm, ym, hm] = measUpdate(xp, yp, hp, sens)
N = size(xp, 2);

pd_w = makedist('Triangular','a',-KC.wbar,'b',1/KC.wbar,'c',KC.wbar);
w = random(pd_w, [4, N]);

beta = fz_xp(xp, yp, hp, sens);

if sum(beta) == 0
    warning('all beta values are zero');
    beta = ones(size(beta))/sum(ones(size(beta)));
else
    beta = beta/sum(beta); %normalize beta
end
[xm, ym, hm] = resample(beta, xp, yp, hp);

end

function [xm, ym, hm] = resample(beta, xp, yp, hp)

N = size(beta, 2);

r = rand(1, N);

for i=1:N
   idx(i) = find(r(i) <= cumsum(beta), 1, 'first');
end

xm = xp(:, idx);
ym = yp(:, idx);
hm = hp(:, idx);

end

function p = fz_xp(xp, yp, hp, z_bar)

mask = z_bar < inf;

d4 = sqrt(xp.^2 + yp.^2);
d3 = sqrt(xp.^2 + (yp - KC.L).^2);
d2 = sqrt((xp - 2 * KC.L).^2 + (yp- KC.L).^2);
d1 = sqrt((xp - 2 * KC.L).^2 + yp.^2);

zm_corr = [d1(1, :);
    d2(1, :);
    d3(2, :);
    d4(2, :)];

zm_incor = [d1(2, :);
    d2(2, :);
    d3(1, :);
    d4(1, :)];

p = fw(z_bar(mask) - zm_corr(mask, :))*KC.sbar + fw(z_bar(mask)-zm_incor(mask, :))*(1-KC.sbar);

end

function prob = fw(w)

prob = zeros(size(w));

mask = w < -KC.wbar | w > KC.wbar;
prob(mask) = 0;

mask = w < KC.wbar & w > 0;
prob(mask) = - w(mask) / KC.wbar^2 + 1/KC.wbar;

mask = w < 0 & w > -KC.wbar;
prob(mask) =  w(mask) / KC.wbar^2 + 1/KC.wbar;

prob = prod(prob, 1);

end

function [xp, yp, hp] = priorUpdate(prevPostParticles, act)
    
    x = prevPostParticles.x;
    y = prevPostParticles.y;
    h = prevPostParticles.h;
    
    N = size(x, 2); % Number of particles
    
    pd_vs = makedist('Triangular','a',-KC.vsbar,'b',1/KC.vsbar,'c',KC.vsbar);
    vs = random(pd_vs, [2, N]);
    pd_v = makedist('Triangular','a',-KC.vsbar,'b',1/KC.vsbar,'c',KC.vsbar);
    v = random(pd_v, [2, N]);
    [xp, yp, hp] = q(x, y, h, act, vs, v);
   
end

function [xk_1, yk_1, hk_1] = q(xm, ym, hm, u, vs, v)

uk = u.*(1 + vs);

dot_x = uk .* cos(hm);
dot_y = uk .* sin(hm);

xk_1 = dot_x*KC.ts + xm;
yk_1 = dot_y*KC.ts + ym;
hk_1 = hm;

%v = makedist('Beta', 'a', 2, 'b', 1);
%KC.vbar
iter = 1;
while(true)
% There was a bounce between the timestep

% Incase there was a bounce with vertical wall during time step
mask1 = (xk_1 < 0 | xk_1 > KC.L * 2);
xk_1(mask1) = - xk_1(mask1) + KC.L * 4 * (xk_1(mask1) > KC.L * 2 ); % without noise
hm(mask1) = pi() - hm(mask1);

% Incase there was a bounce with horizontal wall during time step
mask2 = (yk_1 < 0 | yk_1 > KC.L);
xk_1(mask2) = - xk_1(mask2) + KC.L * 2 * (yk_1(mask2) > KC.L );
hm(mask2) = pi() - hm(mask2);

if sum(mask1 + mask2) == 0
    break;
end
if iter > 100
    warning('maximum iteration exceeded in bounce dynamics');
    break;
iter = iter+1;
end

end
end

function [xm_p, ym_p, hm_p] = roughening(xm, ym, hm)

N = size(xm, 2);
K = 0.06;
d = 6;

E_ix = [max(xm(1, :)) - min(xm(1, :));
        max(xm(2, :)) - min(xm(2, :))];
E_iy = [max(xm(1, :)) - min(xm(1, :));
        max(xm(2, :)) - min(xm(2, :))];
E_ih = [max(xm(1, :)) - min(xm(1, :));
        max(xm(2, :)) - min(xm(2, :))];


sigma_x = K*E_ix*N^(-1/d);
sigma_y = K*E_iy*N^(-1/d);
sigma_h = K*E_ih*N^(-1/d);


randnum = randn(6, N);
delta_x = sigma_x.*randnum(1:2, :);
delta_y = sigma_y.*randnum(3:4, :);
delta_h = sigma_h.*randnum(5:6, :);

xm_p = xm + delta_x;
ym_p = ym + delta_y;
hm_p = hm + delta_h;

iter = 1;
while(true) % Check if roughening results in pushing the robot out of space
    mask_x = xm_p(1, :) > 2* KC.L | xm_p(1, :) < 0 | xm_p(2, :) > 2* KC.L | xm_p(2, :) < 0;
    mask_y = ym_p(1, :) > KC.L | ym_p(1, :) < 0 | ym_p(2, :) > KC.L | ym_p(2, :) < 0;
    
    if sum(sum(mask_x))+sum(sum(mask_x)) ==0
        break;
    else
        
        xm_p(:, mask_x) = xm(:, mask_x) + sigma_x .* randn(2, sum(mask_x));
        ym_p(:, mask_x) = ym(:, mask_x) + sigma_y .* randn(2, sum(mask_x));
    end
    
    if iter > 100
        warning('maximum iterations passed in roughening');
        break;
    end
    iter = iter+1;
    
end


end
