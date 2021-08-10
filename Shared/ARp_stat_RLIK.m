%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARp_stat_RLIK.m                                                         %
%    Calculate likelihood of stationary AR(p) process with restricted     %
%    likelihood                                                           %
% Author: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com> (2019)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val, phi, sigR, Q] = ARp_stat_RLIK(par, y, p)
    % Extract length and number of time series observed
    nObs = size(y, 2);
    nLen = size(y, 1);
    
    % Change to contrasts
    yMean = mean(y, 2);
    for i  = 1:nObs
        y(:,i) = y(:,i)-yMean;
    end
    
    % Local function
    parInv   = @(x)((exp(x)-1)/(exp(x)+1));        
    
    % Extract PACF
        psi = zeros(p,1);
        for i = 1:p
            psi(i) = parInv(par(i));
        end    
    
        % Compute autocorrelation parameters
        allPhi = zeros(p,p);
        for i = 1:p
            allPhi(i,i) = psi(i);
        end
        for k = 2:p
            for j = fliplr(1:(k-1))
                allPhi(k,j) = allPhi(k-1,j) - allPhi(k,k)*allPhi(k-1, k-j);
            end
        end
    
    phi = allPhi(p,:);
    sigM = exp(par(p+1));
    
    % Compute autocorrelations
    if p == 1
        gamma = phi;
    else
        iVec = repmat((1:p)', [1 p]);
        iVec = iVec';
        iVec = iVec(:);
        
        jVec = [1; (1:(p-1))'];
        vVec = [-1; phi(2:p)'];
        for i = 2:p
            jVec = [jVec; flipud((1:i)');(1:(p-i))'];
            phiTmp = phi;
            phiTmp(i) = [];
            vVec = [vVec;-1; phiTmp'];
        end
        Ab = sparse(iVec,jVec,vVec);
        rhs = -phi';
        gamma = Ab\rhs;
    end
    
    % Residual standard deviation
    sigR = sqrt(1-sum(phi'.*gamma))*sigM;
    
    % Calculate precision matrix conditional on y_i, i = 1,...,p
    i = repmat((1:(nLen-p))', [1, 1+p])';
    i = i(:);
    j = i + repmat((0:p)', [nLen-p, 1]);
    v = repmat([-fliplr(phi)';1], [1, nLen-p]);
    v = v(:);
    A = sparse(i,j,v,nLen-p,nLen);
    Q = A'*A/sigR^2;
    
    % Calculate distribution for y_i, i = 1,...,p
    SigStat = toeplitz([1;gamma(1:end-1)])*sigM^2;
    
    % Calculate determinant part of likelihood
    val = -0.5*nLen*(nObs-1)*log(2*pi);
    val = val - 0.5*nLen*log(nObs);
    val = val - 0.5*(nObs-1)*((nLen-p)*log(sigR^2) +log(det(SigStat)));
    %val = -nObs/2*(nLen*log(2*pi) + (nLen-2)*log(sigR^2) - log(det(QStat)));
    
    % Calculate sum of squares part
    for i = 1:nObs
        val = val - 0.5*y(:,i)'*Q*y(:,i) - 0.5*y(1:p,i)'*(SigStat\y(1:p,i));
    end
    
    % Minus likelihood
    val = -val;
    
    if nargout == 4
        Q(1:p,1:p) = Q(1:p, 1:p) + SigStat^-1;
    end
end