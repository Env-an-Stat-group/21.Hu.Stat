%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ProcessData.m                                                           %
%    Read netcdf files, aggregate to annual level and store as .mat       %
%    object.                                                              %
% Author: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com> (2018)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lon, lat, dataV] = ProcessData(runVec, fold, varName)
    % Initialize variable to hold data
    dataV = zeros([320 384 1140 length(runVec)]);

    % Iterate through ensembles
    fNameS = fullfile(fold, 'b.e11.BRCP85C5CNBDRD.f09_g16.');
    fNameE1 = strcat(['.', varName, '.200601-208012.nc']);
    fNameE2 = strcat(['.', varName, '.208101-210012.nc']);
    fNameE3 = strcat(['.', varName, '.200601-210012.nc']);
    for i = runVec
        % Make file name for current ensemble member
        if(i < 34)
            fName1 = strcat([fNameS, sprintf('%03d', i), fNameE1]);
            fName2 = strcat([fNameS, sprintf('%03d', i), fNameE2]);
        elseif(i < 36)
            fName1 = strcat([fNameS, sprintf('%03d', i), fNameE3]);
            fName2 = nan;
        else
            fName1 = strcat([fNameS, sprintf('1%02d', i-35), fNameE3]);
            fName2 = nan;
        end
        
        % Read spatial coordinates
        lat = ncread(fName1, 'TLAT');
        lon = ncread(fName1, 'TLONG');

        % Read data variable
        if(i < 34)
            TREFHT1 = ncread(fName1, 'SST');
            TREFHT2 = ncread(fName2, 'SST');
        else
            TREFHTtmp = ncread(fName1, 'SST');
            TREFHT1 = TREFHTtmp(:,:,1,1:900);
            TREFHT2 = TREFHTtmp(:,:,1,901:end);
        end
        
        dataV(:, :, :, i) = reshape(cat(4, TREFHT1, TREFHT2), [320 384 1140]);
    end

    % Convert to Celsius
    %dataV = dataV - 273.15;

    % Shift data
    %dataV = dataV([146:288, 1:145], :, : ,:);
    %lon = [lon(146:end)-360; lon(1:145)];
end
