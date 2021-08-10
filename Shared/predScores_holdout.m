%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predScores_holdout.m                                                    %
%    Calculate prediction score                                           %
% Author: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com> (2020)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MSE, MAE, BIAS, CRPS] = predScores_holdout(idxP, idxO, allObs, Q, S, Qn)
    % Initialize storage
    MSE = zeros([length(idxP), size(allObs,2)]);
    MAE = MSE;
    CRPS = MSE;
    BIAS = MSE;
    
    % Calculate conditional distribution
    Qc = Q + S(idxO,:)'*Qn(idxO, idxO)*S(idxO,:);
    
    % Make transpose matrices
    ST = S';
    
    % Factorize and permute
    p = amd(Qc);
    Lc = chol(Qc(p,p), 'lower');
    Lct = Lc';
    
    % Calculate prediction variances
    QcInv = SPDE.INLA.sparseQinv(Qc);
    pStd = sqrt(diag(S(idxP,:)*QcInv*S(idxP,:)' + Qn(idxP,idxP)\ones(length(idxP))));
    
    % Calculate across all observations
    for idxObs = 1:size(allObs,2)
        muIy = Lct\(Lc\(ST(p,idxO)*Qn(idxO, idxO)*allObs(idxO, idxObs)));
        muIy(p) = muIy;
        yPred = S(idxP,:)*muIy;
    
        % Stationary scores
        MSE(:, idxObs) = (allObs(idxP,idxObs)-yPred).^2;
        MAE(:, idxObs) = abs(allObs(idxP,idxObs)-yPred);
        zS = (allObs(idxP, idxObs)-yPred)./pStd;
        CRPS(:, idxObs) = -pStd.*(1/sqrt(pi)-2*normpdf(zS)-zS.*(2*normcdf(zS)-1));
        BIAS(:, idxObs) = yPred - allObs(idxP,idxObs);
    end

end