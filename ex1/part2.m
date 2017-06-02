% parameters
N = 50;                     % the number of test
edges = 0:0.1:3;            % bin edge for histogram
Qv = [0.01, 0.1, 1];    

debug = true;

% rms
rms = zeros(N, 3);

for j = 1:size(Qv, 2)
    
    if debug
        disp('===============================================================')
        disp(['Qv = ', num2str(Qv(j))])
    end
    
    % const
    simConst = SimulationConstants();

    % const
    estConst = EstimatorConstants();

    % change Qv
    estConst.VelocityInputPSD = Qv(j);
    
    for i = 1:N
        if debug
            fprintf('(%4d / %4d) \n', i, N);
            tic
        end
        
        % TODO CHECK SEED
        rms(i, 1) = run(simConst, estConst, false, 0);
        
        if debug
            toc;
        end
    end
end

% histogram 
for j = 1:size(Qv, 2)
    fig = histogram(rms(:, j), edges);
    save(fig, ['hist', num2str(j), '.png'])
end