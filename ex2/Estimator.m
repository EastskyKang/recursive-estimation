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
%                       sens(1): distance reported by sensor 1 (metres)
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
N = 500; % obviously, you will need more particles than 10.
if (init)
    % Do the initialization of your estimator here!
    % These particles are the posterior particles at discrete time k = 0
    % which will be fed into your estimator again at k = 1
    
    % start point 
    % A (1st row) : 0 (S1) / 1 (S2)     (same prob)
    % B (2nd row) : 0 (S4) / 1 (S3)     (same prob)
    start_point = rand(2, N) > 0.5;
    
    % convert to coordinates
    x0 = [ones(1, N) * KC.L * 2; zeros(1, N)];
    y0 = start_point * KC.L;
    
    % start orientation      
    h0 = rand(2, N) * (pi()/2); 
    
    % A : [pi/2, pi] (S1) / [-pi, -pi/2] (S2)
    h0(start_point(1,:))  = h0(start_point(1,:)) + (- pi());    % S2
    h0(~start_point(1,:)) = h0(~start_point(1,:)) + (pi()/2);   % S1
    
    % B : [0, pi/2] (S4) / [-pi/2, 0] (S3)
    h0(start_point(2,:))  = h0(start_point(2,:)) + (- pi()/2);  % S3
    h0(~start_point(2,:)) = h0(~start_point(2,:));              % S4
    
    % Replace the following:
    postParticles.x = x0;
    postParticles.y = y0;
    postParticles.h = h0;
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If init = 0, we perform a regular update of the estimator.

% Implement your estimator here!

% -------------------------------------------------------------------------
% S1: prior update
% -------------------------------------------------------------------------

x_m = prevPostParticles.x;
y_m = prevPostParticles.y;
h_m = prevPostParticles.h;

% noise sampling
vs = zeros(2, N);

% f(vs) = 
%   1/vbar_s^2 x + 1/vbar_s     [-vbar_s, 0]
%   - 1/vbar_s^2 x + 1/vbar_s   [0, vbar_s]
%   0                           otherwise
u = rand(2, N);

% case 1 : 0 <= u < 0.5
% vs^2 / (2 * vbar_s^2) + vs / vbar_s + 1/2 - u = 0
% vs = - vbar_s + vbar_s sqrt(2u)
mask = (u < 0.5);
vs(mask) = - KC.vsbar + KC.vsbar * sqrt(u(mask) * 2);  

% case 2 : 0.5 <= u < 1
% - vs^2 / (2 * vbar_s^2) + vs / vbar_s + 1/2 - u = 0
% vs = vbar_s - vbar_s sqrt(2-2u)
mask = (u >= 0.5);
vs(mask) = KC.vsbar - KC.vsbar * sqrt(2 - 2 * u(mask));

% particle propagate
[x_p, y_p, h_p] = q(x_m, y_m, h_m, vs, act);

% -------------------------------------------------------------------------
% S2: measurement update
% -------------------------------------------------------------------------

% distance
d1_A = sqrt((x_p(1,:) - KC.L * 2).^2 + (y_p(1,:)).^2);          % (1 x N)
d1_B = sqrt((x_p(2,:) - KC.L * 2).^2 + (y_p(2,:)).^2);          % (1 x N)
d2_A = sqrt((x_p(1,:) - KC.L * 2).^2 + (y_p(1,:) - KC.L).^2);   % (1 x N)
d2_B = sqrt((x_p(2,:) - KC.L * 2).^2 + (y_p(2,:) - KC.L).^2);   % (1 x N)
d3_A = sqrt((x_p(1,:)).^2 + (y_p(1,:) - KC.L).^2);              % (1 x N)
d3_B = sqrt((x_p(2,:)).^2 + (y_p(2,:) - KC.L).^2);              % (1 x N)
d4_A = sqrt((x_p(1,:)).^2 + (y_p(1,:)).^2);                     % (1 x N)
d4_B = sqrt((x_p(2,:)).^2 + (y_p(2,:)).^2);                     % (1 x N)

d_correct = [d1_A; d2_A; d3_B; d4_B];
d_wrong = [d1_B; d2_B; d3_A; d4_A];

% measurement
z = sens;      % (4 x 1)      

% remove invalid measurements 
mask = ~isinf(z);

z = z(mask, :);
d_wrong = d_wrong(mask, :);
d_correct = d_correct(mask, :);

% calculate f(z|p) = f_w(z - d_B) sbar + f_w(z - d_A) (1 - sbar)
beta = fw(z - d_wrong) * KC.sbar ...
    + fw(z - d_correct) * (1 - KC.sbar);        % (a x N)

% columnwise multiplication (w1, w2, w3, w4 are indep)
beta = prod(beta, 1);                           % (1 x N)

if sum(beta) == 0
    % sum of beta is zero (beta is all zero)
    % measurement model cannot deal with large value of measurement error
   
    warning('beta is all zero!')
    beta = ones(1, N) * 1 / N;                  % TODO is there better way?
end

% normalize beta
beta = beta / sum(beta);                        %(1 x N)

% -------------------------------------------------------------------------
% resampling

if ~isempty(z)
    % accumulated beta
    cum_beta = cumsum(beta);
    
    % random numbers
    u = rand(1, N);
    
    % idx of particles
    idx = zeros(1, N);
    
    for i = 1:N
        idx(i) = find(u(i) <= cum_beta, 1);
    end
    
    x_m = x_p(:, idx);
    y_m = y_p(:, idx);
    h_m = h_p(:, idx);

else
    % if no measurement (all elm of sens are inf)
    x_m = x_p;
    y_m = y_p;
    h_m = h_p;
end

% -------------------------------------------------------------------------
% roughening
[x_m, y_m, h_m] = Roughening(x_m, y_m, h_m, N);

% Replace the following:
postParticles.x = x_m;
postParticles.y = y_m;
postParticles.h = h_m;

end % end estimator


function [x, y, h] = q(x, y, h, v, u)
   
    % x     x_m[k-1]    (2 x N)
    % y     y_m[k-1]    (2 x N)
    % h     h_m[k-1]    (2 x N)
    % v     v_s[k-1]    (2 x N)
    % u     u[k-1]      (2 x 1)
    
    % debug
    debug = false;
    
    if (debug)
        figure(10)
        plot(x(1,:), y(1,:), 'b.');
        hold on 
        plot(x(2,:), y(2,:), 'r.');
        hold off
        drawnow
        pause(1)
    end
    
    % the number of particle
    n = size(x, 2);
        
    % u(t)
    u = u .* (1 + v);
    
    % time remains for each particles 
    time_remains = ones(2, n) * KC.ts;
    
    while nnz(time_remains) ~= 0
        % update position until every remain time is zero
        % this code vectorize inputs for fast run. 
        
        % velocity (column = particle)
        xdot = u .* cos(h);      % 2 x N
        ydot = u .* sin(h);      % 2 x N
        
        % -----------------------------------------------------------------
        % time to bounce (calculate to every wall)
        
        % mask (bounce can happens to that wall ?)
        % if mask = 0 then time to wall is inf (or NaN)
        mask = (xdot < 0);
        time_to_lt = (0 - x) ./ (xdot .* mask);         % to left wall   (2 x N)
        time_to_lt = abs(time_to_lt);
        
        mask = (xdot > 0);
        time_to_rt = (KC.L * 2 - x) ./ (xdot .* mask);  % to right wall  (2 x N)
        time_to_rt = abs(time_to_rt);
        
        mask = (ydot > 0);
        time_to_up = (KC.L * 1 - y) ./ (ydot .* mask);  % to upper wall  (2 x N)
        time_to_up = abs(time_to_up);
        
        mask = (ydot < 0);
        time_to_lw = (0 - y) ./ (ydot .* mask);         % to lower wall  (2 x N)
        time_to_lw = abs(time_to_lw);
        
        % now time_to_x contains time to bounce 
        % if its value less than time_remains then bounce will happens
        
        % matrix time to bounce contains followings (2 x N)
        %
        % case 1 : if bounce will happens then minimum time to bounce 
        %          (wall_idx = 1, 2, 3, 4)
        % case 2 : if bounce will not happens then time remains
        %          (wall_idx = 5)
        time_to_bounce = cat(3, time_to_lt, time_to_rt, time_to_up, time_to_lw, time_remains);
        [time_to_bounce, wall_idx] = min(time_to_bounce, [], 3);

        % -----------------------------------------------------------------
        % update x and y  
        % if time_remains is zero, it does not update
        x = x + xdot .* time_to_bounce;   
        y = y + ydot .* time_to_bounce;
        
        % -----------------------------------------------------------------
        % update h
        % it only for wall_idx = 1, 2, 3, 4
        
        % generate random numbers 
        % f(v) = cv^2 for [-v_bar, v_bar]
        % c is constant
        c = 3 / (2 * KC.vbar)^3;
        
        % by F(v) = c/3 v^3 + c/3 v_bar^3
        % u is random number [0, 1]
        u = rand(2, n);
        
        % v = (3u / c - v_bar^3)^1/3
        v = (3 * u / c - KC.vbar^3).^ 1/3;       % (2 x N)
        
        % ideal post-bound angle (to every wall) 
        post_h_lt = atan2(ydot, -xdot);                 
        post_h_rt = atan2(ydot, -xdot);
        post_h_up = atan2(-ydot, xdot);
        post_h_lw = atan2(-ydot, xdot);
        
        % actual post-bound angle (to every wall)
        post_h_lt = post_h_lt .* (1 + v);
        post_h_rt = post_h_rt .* (1 + v);
        post_h_up = post_h_up .* (1 + v);
        post_h_lw = post_h_lw .* (1 + v);
        
        h(wall_idx == 1) = post_h_lt(wall_idx == 1);
        h(wall_idx == 2) = post_h_rt(wall_idx == 2);
        h(wall_idx == 3) = post_h_up(wall_idx == 3);
        h(wall_idx == 4) = post_h_lw(wall_idx == 4);
        
        % -----------------------------------------------------------------
        % update time remains
        time_remains = time_remains - time_to_bounce;
        
        if (debug)
            figure(10)
            plot(x(1,:), y(1,:), 'b.');
            hold on
            plot(x(2,:), y(2,:), 'r.');
            hold off
            drawnow
            pause(1)
        end
    end
end

function [p] = fw(w)
    % w     (w1[k], w2[k], w3[k], w4[k])'   (4 x N)
   
    % the number of particles
    p = zeros(size(w));
    
    if size(w, 1) == 0
        return
    end
    
    % case 1 : w >= w_bar
    mask = (KC.wbar < w);
    p(mask) = 0;
    
    % case 2 : w_bar > w >= 0
    mask = (0 <= w & w < KC.wbar);
    p(mask) = - w(mask) / KC.wbar^2 + 1/ KC.wbar;
    
    % case 3 : 0 > w >= - w_bar
    mask = (-KC.wbar <= w & w < 0);
    p(mask) = w(mask) / KC.wbar^2 + 1/ KC.wbar;
    
    % case 4 : - w_bar > w
    mask = (w < -KC.wbar);
    p(mask) = 0;
    
end

function [x_m, y_m, h_m] = Roughening(x_m, y_m, h_m, N)
        
    K = 0.04;        % tuning parameter
    d = 6;          % dimension of state space
    
    % inter-sample variability
    Ex_A = max(x_m(1,:)) - min(x_m(1,:));
    Ex_B = max(x_m(2,:)) - min(x_m(2,:));
    Ey_A = max(y_m(1,:)) - min(y_m(1,:));
    Ey_B = max(y_m(2,:)) - min(y_m(2,:));
    Eh_A = max(h_m(1,:)) - min(h_m(1,:));
    Eh_B = max(h_m(2,:)) - min(h_m(2,:));
    
    sig_x = K * [Ex_A; Ex_B] * (N ^ (-1/d));
    sig_y = K * [Ey_A; Ey_B] * (N ^ (-1/d));
    sig_h = K * [Eh_A; Eh_B] * (N ^ (-1/d));
    
    % sampling mask for particles
    % if mask is true then resample
    mask_x = true(size(x_m));
    mask_y = true(size(y_m));
    
    % delta
    delta_x = zeros(size(x_m));
    delta_y = zeros(size(y_m));
    delta_h = normrnd(zeros(size(h_m)), sig_h .* ones(size(h_m)));
    
    while ~((nnz(mask_x) == 0) && (nnz(mask_y) == 0))
        % until valid... (particles never goes over the walls)
        
        % sample from normal distribution
        delta_x = randn(size(x_m)) .* sig_x .* mask_x + delta_x .* (~mask_x);
        delta_y = randn(size(y_m)) .* sig_y .* mask_y + delta_y .* (~mask_y);
       
        % update mask 
        mask_x = ~((0 <= (x_m + delta_x)) & ((x_m + delta_x) <= (KC.L * 2)));
        mask_y = ~((0 <= (y_m + delta_y)) & ((y_m + delta_y) <= (KC.L * 1)));
        
    end
    
    x_m = x_m + delta_x;
    y_m = y_m + delta_y;
    h_m = h_m + delta_h;
end