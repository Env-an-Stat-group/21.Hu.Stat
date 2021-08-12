%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference.m                                                             %
%    Read and prepare data, and fit model, .                              %
%    Use restricted likelihood approach for time and likelihood for space %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setenv('LD_LIBRARY_PATH', '../Shared/:/opt/crc/g/gsl/2.5/gcc/lib/:/afs/crc.nd.edu/x86_64_linux/r/R/3.6.2/gcc/4.8.5/lib64/R/lib/');

%% Setup
    % Number of workers
    parpool(24);
    
    % General code
    addpath('../Shared')
  
%% Load sample data
    % Create datasubset for training stochastic generator
        dim1 = 320;
        dim2 = 384;
        disp('Loading mesh:')
        load('../../Data/mesh.mat');
        disp('Loading observations:')
        load('../../Data/allObs.mat');
        
    %% Stationary model
        %% Fit all data
             disp('Fitting stationary model. All data:');
             sTime = tic;
             
             % Create SPDE model
             tWindow = 1:size(allObs,2);
             tensorParOld = 1;
             OptStat = SPDE.Optimizer.makeStationaryModel(vLoc, tt, tv, loc, allObs(:, tWindow), 0, tensorParOld);

             % Load starting value
             x0 = [-5.6; 0; 0; -1; 3.51;2.68];

             % Fit model
             fitStat = true;
             if(fitStat)
                % Set optimization function
                fun = @(par)(OptStat.logLikelihood(par, [], 1e-4, [], 1, sqrt(eps), 0));

                %% Optimize
                [xStat, valS] = fminunc(fun, x0, optimset('MaxIter', 2000, 'Display', 'iter-detailed', 'GradObj', 'on', 'LargeScale', 'off'));

                % Store result
                tFitStat = toc(sTime);
                save('../../Results/spatial_stat_single.mat', 'xStat', 'tFitStat');
            else
                load('../../Results/spatial_stat_single.mat');
            end

            toc;
            
                        
     %% Non-stationary model
        %% All data
            disp('Fitting non-stationary model. All data')
            sTime = tic;

            % Create time window
            %tWindow = 1:(10*12*5);
            tWindow = 1:size(allObs,2);
            
            % Create SPDE model
            OptNStat = SPDE.Optimizer.makeNonStatModel(vLoc, tt, tv, loc, allObs(:, tWindow), 0, 1);

            % Start value
            x0 = rand(78, 1)*0.1-0.05;
            x0(1) = xStat(1)+xStat(4)/0.28;
            x0(2) = xStat(1);
            x0(end-1) = xStat(end-1);
            x0(end) = xStat(end);
            
            % Iterate through windows
            fitNStat = true;
            if(fitNStat)
                % Set OptStatimization function
                fun = @(par)(OptNStat.logLikelihood(par, [], 1e-4, [], 1, sqrt(eps), 0));

                %% Optimize
                [xNStat, valNS] = fminunc(fun, x0, optimset('MaxIter', 600, 'Display', 'iter-detailed', 'GradObj', 'on', 'LargeScale', 'off', 'TolFun', 1e-5));

                tFitNStat = toc(sTime);
                save('../../Results/spatial_nstat_single.mat', 'xNStat', 'tFitNStat');
            else
                load('../../Results/spatial_nstat_single.mat');
            end
            toc(sTime);
            
    %% Fit wind-based models
        %% Three extra parameters
            disp('Fitting wind model. All data:');
            sTime = tic;

            % Create SPDE model
            tWindow = 1:size(allObs,2);
            tensorParOld = 1;
            OptWStat = SPDE.Optimizer.makeStationaryModel(vLoc, tt, tv, loc, allObs(:, tWindow), 0, tensorParOld);
            OptWStat.addWind(1);

            % Load starting value
            x0 = rand(7, 1)*0.1-0.05;

            % Fit model
            fitWStat = true;
            if(fitWStat)
                % Set optimization function
                fun = @(par)(OptWStat.logLikelihood(par, [], 1e-4, [], 1, sqrt(eps), 0));

                %% Optimize
                [xWStat3, valWS3] = fminunc(fun, x0, optimset('MaxIter', 2000, 'Display', 'iter-detailed', 'GradObj', 'on', 'LargeScale', 'off'));

                % Store result
                tFitWStat3 = toc(sTime);
                save('../../Results/spatial_Wstat3_single.mat', 'xWStat3', 'tFitWStat3');
            else
                load('../../Results/spatial_Wstat3_single.mat');
            end

            toc;
            
        %% Spatially varying range scaling 3 par
            disp('Fitting wind model; more complex. All data:');
            sTime = tic;

            % Create SPDE model
            tWindow = 1:size(allObs,2);
            tensorParOld = 1;
            OptWNStat = SPDE.Optimizer.makeNonStatModel(vLoc, tt, tv, loc, allObs(:, tWindow), 0, tensorParOld, [4 0 4]);
            OptWNStat.addWind(1);

            % Load starting value
            x0 = rand(81, 1)*0.1-0.05;
            x0(1:76) = xNStat(1:76);
            x0(80:81) = xNStat(77:78);
            
            % Fit model
            fitWNStat = true;
            if(fitWNStat)
                % Set optimization function
                fun = @(par)(OptWNStat.logLikelihood(par, [], 1e-4, [], 1, sqrt(eps), 0));

                %% Optimize
                [xWNStat3, valWNS3] = fminunc(fun, x0, optimset('MaxIter', 400, 'Display', 'iter-detailed', 'GradObj', 'on', 'LargeScale', 'off'));

                % Store result
                tFitWNStat3 = toc(sTime);
                save('../../Results/spatial_WNstat3_single.mat', 'xWNStat3', 'tFitWNStat3');
            else
                load('../../Results/spatial_WNstat3_single.mat');
            end

            toc;

        %% Add neural net wind
            %% Spatially varying range harmonics and Neural Network + harmonics
                disp('Fitting NN model');
                disp('Load estimators from Non-stat model as starting points')
                load('../../Results/spatial_nstat_single.mat','xNStat')
                sTime = tic;

                % Create SPDE model
                tWindow = 1:size(allObs,2);
                tensorParOld = 1;
                OptNN = SPDE.Optimizer.makeNN(vLoc, tt, tv, loc, allObs(:, tWindow), 0, tensorParOld, [4 0 4]);
                nbNum=4;
                nodeNum=3;
                OptNN.addWindNN(nbNum, nodeNum);
                % Load starting value
                size = 76 + ((nbNum + 1) * nodeNum + 2 * nodeNum + 1)*3 + 2;
                x0 = rand(size, 1)*0.1-0.05;
                x0(1:76) = xNStat(1:76);
                x0((size-1):size) = xNStat(77:78);

                % Iterate through windows
                fitNN = true;
                if(fitNN)
                    fun = @(par)(OptNN.logLikelihood(par, [], 1e-4, [], 1, sqrt(eps), 0));

                    %% Optimize
                    [xNN, valNN] = fminunc(fun, x0, optimset('MaxIter', 200, 'Display', 'iter-detailed', 'GradObj', 'off', 'LargeScale', 'off'));

                    % Store result
                    tFitWNStat2 = toc(sTime);
                    save('../../Results/spatial_NNhm.mat', 'xNN','valNN','tFitWNStat2');
                else
                    load('../../Results/spatial_NNhm.mat', 'xNN','valNN','tFitWNStat2');
                end
                toc(sTime)
