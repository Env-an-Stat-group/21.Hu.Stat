%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FitVARpprofile_RLIK.m                                                   %
%    Fit VAR2 process separately at each pixel with restricted likelihood %
% Author: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com> (2019)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi, Sig, parMat, optRes] = FitVARpprofile_RLIK(dataO, p)
    % Initialize problem
    d1 = size(dataO, 1);
    d2 = size(dataO, 2);
    Phi = zeros(d1, d2, p);
    Sig  = zeros(d1, d2);
    optRes = Sig*nan;
    parMat = zeros([d1 d2 p+1]);

    % Iterate through each pixel
    tmpSize = size(dataO, 4);
    parfor i = 1:d1
        for j = 1:d2
            if(sum(isnan(dataO(i,j,:))) > 0)
                Phi(i,j,:) = NaN;
                Sig(i,j) = NaN;
                parMat(i,j,:) = NaN;
                continue;
            end
            % Preprocess mean structure
            yTime = reshape(dataO(i, j, :), [1140, tmpSize]);
           
            % Estimate AR(p) structure
            x0 = zeros([p+1 1]);
            optFun = @(x)(ARp_stat_RLIK(x, yTime, p));
            [x, ~, optRes(i, j)] = fminunc(optFun, x0, optimset('Display', 'off'));
            [~, Phi(i,j,:), Sig(i,j)] = optFun(x);
            parMat(i,j,:) = x;
        end
    end
end