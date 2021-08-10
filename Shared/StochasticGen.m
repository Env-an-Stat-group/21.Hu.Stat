%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% StochasticGen.m                                                         %
%    Generate new realizations from stochastic generator                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ensemble] = StochasticGen(muAll, Phi, Sig, pOrder, S, Q, Nrun, Nburn, seedVal, iid)
    % Set seed
    
    rng(seedVal);
    
    % Precompute
    p = amd(Q);
    Lt = chol(Q(p,p), 'lower')';

    % Create matrices describing VAR(pOrder)
    nDim = 86212;
    Dphi = cell(pOrder,1);
    for i = 1:pOrder
        temp_Phi = squeeze(Phi(:,:,i));
        temp_Phi = temp_Phi(:);
        temp_Phi = temp_Phi(~isnan(temp_Phi));
        Dphi{i} = spdiags(temp_Phi, 0, nDim, nDim);
    end
    temp_Sig = Sig(:);
    temp_Sig = temp_Sig(~isnan(temp_Sig));
    Dsig  = spdiags(temp_Sig, 0, nDim, nDim);

    % Initialize storage
    % Initialize storage for monthly average for two location
    refLoc = [35007, 29159]; % Atlantic; Indian
    if Nrun > 1
        ensemble = zeros([2, 1140, Nrun]);
    else
        ensemble = zeros([size(S,1), 1140, Nrun]);
    end

    % Generate runs
    Ntot = 1140+Nburn;
    for iRun = 1:Nrun
        disp('Generating run number');
        disp(iRun);

        % Simulate innovations
        if(iid == 0)
            eMat = Lt\randn(size(Q,1), Ntot);
            eMat(p,:) = eMat;

            % Project to grid
            eMat = S*eMat;
        else
            eMat = randn([size(S,1), Ntot]);
        end

        % Generate VAR(2) process
        eMat = Dsig*eMat;
        for it = (pOrder+1):Ntot
            for pLag = 1:pOrder
                eMat(:,it) = Dphi{pLag}*eMat(:,it-pLag) + eMat(:,it);
            end
        end
 
        % Extract one run
        muAll(muAll == 0) = NaN;
        temp_mu = reshape(muAll,[320*384, size(muAll, 3)]);
        nanvec = sum(isnan(temp_mu), 2);
        isobs = nanvec == 0;
        temp_mu = temp_mu(isobs,:);
        if Nrun > 1
            ensemble(1,:,iRun) = temp_mu(refLoc(1),:) +eMat(refLoc(1),(Ntot-1139):Ntot);
            ensemble(2,:,iRun) = temp_mu(refLoc(2),:) +eMat(refLoc(2),(Ntot-1139):Ntot);
        else
            ensemble(:,:,iRun) = temp_mu +eMat(:,(Ntot-1139):Ntot);
        end
        
    end
end