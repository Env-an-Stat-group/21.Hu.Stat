%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GetControlRun.m                                                         %
%    Read control run                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lon, lat, dataV, dataJan] = GetControlRun(refS, refE, fold)
    % Initialize variable to hold data
    dataV = zeros([288 192 100]);

    % File name
    fName = fullfile(fold, 'b.e11.B1850C5CN.f09_g16.005.cam.h0.TREFHT.200001-209912.nc');
    
    % Read spatial coordinates
    lat = ncread(fName, 'lat');
    lon = ncread(fName, 'lon');

    % Read time
    time = ncread(fName, 'time');

    % Length of months
    mLen = diff([0; time(1:12)]);
    mLen(1) = 31;

    % Read data variable and aggregate to annual mean
    TREFHT = ncread(fName, 'TREFHT');

    for yr = 1:100
        for mth = 1:12
            cVal = TREFHT(:,:,12*(yr-1)+mth);
            dataV(:,:,yr) = dataV(:,:,yr) + mLen(mth)*cVal;
        end
        dataV(:, :, yr) = dataV(:, :, yr)/sum(mLen);
    end

    % Convert to Celsius
    dataV = dataV - 273.15;

    % Shift data
    dataV = dataV([146:288, 1:145], :, (refS-1999):(refE-1999));
    lon = [lon(146:end)-360; lon(1:145)];
    
    % Get January average
    dataJan = TREFHT([146:288, 1:145], :, 1 + 12*((refS-1999):(refE-1999)))-273.15;
end

