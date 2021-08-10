%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARp_RLIK.m                                                              %
%    Simple AR(p) estimation using multiple realizations and restricted   %
%    likelihood.                                                          %
% Author: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com> (2019)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val, phi, sigR, llik] = ARp_RLIK(par, y, order)
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
    psi = zeros(12,1);
    for i = 1:order
        psi(i) = parInv(par(i));
    end    
    sigR = exp(par(13));
    
    % Compute autocorrelation parameters
    allPhi = zeros(12,12);
    for i = 1:12
        allPhi(i,i) = psi(i);
    end
    for k = 2:12
        for j = fliplr(1:(k-1))
            allPhi(k,j) = allPhi(k-1,j) - allPhi(k,k)*allPhi(k-1, k-j);
        end
    end
    
    phi = allPhi(order,:);
    
    % Calculate precision matrix
    i = reshape(repmat(1:(nLen-12), [13 1]), [(nLen-12)*13, 1]);
    j = reshape([1:(nLen-12);
                 2:(nLen-11);
                 3:(nLen-10);
                 4:(nLen-9);
                 5:(nLen-8);
                 6:(nLen-7);
                 7:(nLen-6);
                 8:(nLen-5);
                 9:(nLen-4);
                 10:(nLen-3);
                 11:(nLen-2);
                 12:(nLen-1);
                 13:nLen], [(nLen-12)*13 1]);
    v = reshape(repmat([-fliplr(phi), 1], [nLen-12, 1])', [(nLen-12)*13 1]);
    A = sparse(i, j, v);
    Q = A'*A;
        
    % Calculate sum of squares
    yTilde = y;
    ssq = 0;
    for i = 1:nObs
        yTilde(:,i) = y(:,i);
        ssq = ssq + yTilde(:,i)'*Q*yTilde(:,i);
    end
    likeSSQ = -0.5*ssq/sigR^2;

    % Calculate determinant part
    likeDet = -0.5*(nLen-12)*((nObs-1)*(log(2*pi) + log(sigR^2))+log(nObs));
    
    % Minus likelihood
    llik = likeSSQ+likeDet;
    val = -llik;
end