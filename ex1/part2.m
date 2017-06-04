% parameters
N = 50;                     % the number of test
edges = 0:0.1:3;            % bin edge for histogram

Qv_sim = [0.01, 0.1, 1];
Qv_est = [0.01, 0.1, 1];    

debug = true;
save_out = true;

%% RMS
rms = zeros(N, 3);

for j = 1:size(Qv_est, 2)
    
    if debug
        disp('===============================================================')
        disp(['Qv = ', num2str(Qv_est(j))])
    end
    
    % const
    simConst = SimulationConstants();

    % const
    estConst = EstimatorConstants();

    % change Qv
    estConst.VelocityInputPSD = Qv_sim(j);
    estConst.VelocityInputPSD = Qv_est(j);
    
    for i = 1:N
        if debug
            fprintf('(%4d / %4d) \n', i, N);
            tic
        end
        
        rms(i, j) = run(simConst, estConst, false, 0);
        
        if debug
            toc;
        end
    end
end

%% RMS HISTOGRAM
disp('===============================================================')
disp('draw and save histogram')

for j = 1:size(Qv_est, 2)
    figure(j)
    fig = histogram(rms(:, j), edges);
    
    title(['RMS error histogram with Qv = ', num2str(Qv_est(j))])
    xlabel('RMS error (Root Mean Squared Error) (m)')
    
    if save_out
        saveas(fig, ['hist', num2str(j), '.png'])
    end
end

%% MEAN AND VARIANCE OF RMS mean and variance 
disp('===============================================================')
disp('calculating mean and variance of RMS')

mean_rms = mean(rms);
var_rms = var(rms);

disp('mean = ')
disp(mean_rms)

if save_out
    save('mean_rms.txt', 'mean_rms', '-ascii', '-tabs', '-double')
end

disp('var = ')
disp(var_rms)

if save_out
    save('var_rms.txt', 'var_rms', '-ascii', '-tabs', '-double')
end