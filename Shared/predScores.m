%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predScores.m                                                            %
%    Calculate prediction scores for stationary and non-stationary        %
% Author: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com> (2019)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MSE, MAE, BIAS, CRPS] = predScores(idxP, idxO, allObs, Qs, Ss, QnS, Qns, Sns, QnNS)
    % Initialize storage
    MSE = zeros([length(idxP), size(allObs,2), 2]);
    MAE = MSE;
    CRPS = MSE;
    BIAS = MSE;
    
    % Calculate conditional distribution
    QcS = Qs + Ss(idxO,:)'*QnS(idxO, idxO)*Ss(idxO,:);
    QcNS = Qns + Sns(idxO,:)'*QnNS(idxO, idxO)*Sns(idxO,:);
    
    % Make transpose matrices
    SsT = Ss';
    SnsT = Sns';
    
    % Factorize and permute
    p = amd(QcS);
    LcS = chol(QcS(p,p), 'lower');
    LcSt = LcS';
    LcNS = chol(QcNS(p,p), 'lower');
    LcNSt = LcNS';

    % Calculate prediction variances
    QcInvS = SPDE.INLA.sparseQinv(QcS);
    QcInvNS = SPDE.INLA.sparseQinv(QcNS);
    pStdS = sqrt(diag(Ss(idxP,:)*QcInvS*Ss(idxP,:)' + QnS(idxP,idxP)\ones(length(idxP))));
    pStdNS = sqrt(diag(Sns(idxP,:)*QcInvNS*Sns(idxP,:)' + QnNS(idxP,idxP)\ones(length(idxP))));
    
    % Calculate across all observations
    for idxObs = 1:size(allObs,2)
        muIyS = LcSt\(LcS\(SsT(p,idxO)*QnS(idxO, idxO)*allObs(idxO, idxObs)));
        muIyS(p) = muIyS;
        yPredS = Ss(idxP,:)*muIyS;
        muIyNS = LcNSt\(LcNS\(SnsT(p,idxO)*QnNS(idxO, idxO)*allObs(idxO, idxObs)));
        muIyNS(p) = muIyNS;
        yPredNS = Sns(idxP,:)*muIyNS;
    
        % Stationary scores
        MSE(:, idxObs, 1) = (allObs(idxP,idxObs)-yPredS).^2;
        MAE(:, idxObs, 1) = abs(allObs(idxP,idxObs)-yPredS);
        zS = (allObs(idxP, idxObs)-yPredS)./pStdS;
        CRPS(:, idxObs, 1) = -pStdS.*(1/sqrt(pi)-2*normpdf(zS)-zS.*(2*normcdf(zS)-1));
        BIAS(:, idxObs, 1) = yPredS - allObs(idxP,idxObs);
        
        % Non-stationary scores
        MSE(:, idxObs, 2) = (allObs(idxP,idxObs)-yPredNS).^2;
        MAE(:, idxObs, 2) = abs(allObs(idxP,idxObs)-yPredNS);
        zNS = (allObs(idxP, idxObs)-yPredNS)./pStdNS;
        CRPS(:, idxObs, 2) = -pStdNS.*(1/sqrt(pi)-2*normpdf(zNS)-zNS.*(2*normcdf(zNS)-1));
        BIAS(:, idxObs, 2) = yPredNS - allObs(idxP,idxObs);
    end

end