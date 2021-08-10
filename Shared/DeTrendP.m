%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DeTrendP.m                                                              %
%    Remove AR(p) structure.                                              %
% Author: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com> (2019)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataD] = DeTrendP(dataO, mu, Phi, Sig, p)
    % Remove mean
    for i = 1:size(dataO, 4)
        dataO(:,:,:,i) = dataO(:,:,:,i) - mu;
    end
    
    % Collapse two first dimensions
    oSize = size(dataO);
    tSize = oSize;
    tSize = [prod(tSize(1:2)), tSize(3:end)]; 
    dataO = reshape(dataO, tSize);
    
    % Create matrices describing VAR(2)
    nDim = tSize(1);
    Dphi = cell(6,1);
    for i = 1:6
        Dphi{i} = spdiags(reshape(Phi(:,:,i), [nDim 1]), 0, nDim, nDim);
    end
    Dsig  = spdiags(Sig(:), 0, nDim, nDim);
    
    % De-correlate
    dataD = dataO(:, (p+1):end, :);
    for eIdx = 1:oSize(4)
        for ord = 1:p
            dataD(:, :, eIdx) = dataD(:, :, eIdx) - Dphi{ord}*dataO(:, (p+1-ord):(end-ord), eIdx);
        end
        dataD(:, :, eIdx) = Dsig\dataD(:, :, eIdx);
    end
end
