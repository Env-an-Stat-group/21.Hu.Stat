classdef Optimizer < handle
	% OPTIMIZER
	%    Calculate log-likelihood and optimize

	properties
        QMaker;
		A;
		obs;
        dMat;
        lMat;
        lMatE;
	end

	methods
		function obj = Optimizer(vLoc, tt, tv, tensorPar)
			% OPTIMIZER
			%    Constructor
            
            if(nargin < 4)
                tensorPar = 1;
            end

            obj.QMaker = SPDE.PrecisionMaker(vLoc, tt, tv, tensorPar);
            
            %% Add distances from border
                % Get distances from continent borders
                S = shaperead('../Data/continents/continent');
                X1 = S(1).X;
                Y1 = S(1).Y;
                for i = 2:7
                    [xa, ya] = polybool('union', X1, Y1, S(i).X, S(i).Y);
                    X1 = xa;
                    Y1 = ya;
                end
                X1 = [X1,S(8).X];
                Y1 = [Y1,S(8).Y];
                
                % Remove line at latitude -90
                idx = (Y1<-90+1e-3);
                X1(idx) = [];
                Y1(idx) = [];

                cBorder = cell(1);
                for i = 1:1
                    [x, y, z] = sph2cart(X1*pi/180, Y1*pi/180, 1);
                    cBorder{i} = [x', y', z'];
                end

                % Get coordinates of centres of vertices
                xCenter = obj.QMaker.meshCLoc;

                % Find distances
                dGrid = abs(inf*xCenter(:,1));
                tic;
                for j = 1:length(dGrid)
                    % Subsample
                    idx = floor(linspace(1, size(cBorder{1}, 1), 10000));
                    dGrid(j) = min(dGrid(j), min(pdist2([xCenter(j,1), xCenter(j,2), xCenter(j,3)], cBorder{1}(idx,:))));
                end
                toc;

                obj.dMat = dGrid;
                
             %% Add land/sea
                 % Find az, el of centres
                [az, el, ~] = cart2sph(obj.QMaker.meshCLoc(:,1), obj.QMaker.meshCLoc(:,2), obj.QMaker.meshCLoc(:, 3));
                lon = az*360/(2*pi);
                lat = el*360/(2*pi);
                
                xv = lon*0;
                for i = 1:8
                    xv = xv + inpolygon(lon, lat, S(i).X, S(i).Y);
                end

                % Interpolate with nearest neighbour
                obj.lMat = xv;
                
                gridLandSeaE = repmat(obj.lMat, 1, 3)';
                gridLandSeaE = gridLandSeaE(:);
                obj.lMatE = gridLandSeaE;
                
                toc;
        end
        
        function addLandSeaCov(obj)
           % ADDLANDSEACOV
           %    Add land sea covariate
           
            % Set basis
            obj.QMaker.basisRho = [obj.QMaker.basisRho, obj.lMat];
            gridLandSeaE = repmat(obj.lMat, 1, 3)';
            gridLandSeaE = gridLandSeaE(:);
            obj.lMatE = gridLandSeaE;
              
            % Set indicies
            currMax = max([obj.QMaker.idxRho, obj.QMaker.idxSigma, obj.QMaker.idxDVx, obj.QMaker.idxDVy, obj.QMaker.idxDVz]);
            obj.QMaker.idxRho = [obj.QMaker.idxRho, currMax+1];
        end
        
        function addWindNN(obj, nbNum, nodeNum)
           % ADDWIND
           %    Add wind covariate
           
            % Load zonal wind
            load('../Data/UBOT.mat');

            % Load meridional wind
            load('../Data/VBOT.mat');

            % Expand longitude and latitude to a mesh
            [gLat, gLon] = meshgrid(lat, lon);
            % Fix edge
            gLon = [gLon;gLon(1,:)+360];
            gLat = [gLat;gLat(1,:)];
            UBOT = [UBOT;UBOT(1,:)];
            VBOT = [VBOT;VBOT(1,:)];
            
            % Re-grid
            [az, el, ~] = cart2sph(obj.QMaker.meshELoc(:,1), obj.QMaker.meshELoc(:,2), obj.QMaker.meshELoc(:, 3));
            lonE = az*360/(2*pi);
            lonE(lonE < 0) = lonE(lonE < 0) + 360;
            latE = el*360/(2*pi);
            Uq = interp2(gLat, gLon, UBOT, latE, lonE);
            Vq = interp2(gLat, gLon, VBOT, latE, lonE);
            
            % Map to cartesian vector field
            cartVec = [cos(el).*(-sin(az)).*Uq-sin(el).*cos(az).*Vq, cos(el).*cos(az).*Uq-sin(el).*sin(el).*Vq, cos(el).*Vq]*pi/180;
           
            % find which grid the edge centers are located at
            dLon = lon(2) - lon(1);
            dLat = lat(2) - lat(1);
            indLon = floor(lonE / dLon) + 1;
            indLat = floor((latE + 90) / dLat) + 1;
            indLon(indLon > length(lon)) = length(lon);
            indLat(indLat > length(lat)) = length(lat);
            % find neiboring grid
            nbDist = 1;
            nbLat1 = indLat - nbDist;
            nbLat1(nbLat1 > length(lat)) = nbLat1(nbLat1 > length(lat)) - length(lat);
            nbLat1(nbLat1 < 1) = nbLat1(nbLat1 < 1) + length(lat);
            nbLat2 = indLat + nbDist + 1;
            nbLat2(nbLat2 > length(lat)) = nbLat2(nbLat2 > length(lat)) - length(lat);
            nbLat2(nbLat2 < 1) = nbLat2(nbLat2 < 1) + length(lat);
            nbLon1 = indLon - nbDist;
            nbLon1(nbLon1 > length(lon)) = nbLon1(nbLon1 > length(lon)) - length(lon);
            nbLon1(nbLon1 < 1) = nbLon1(nbLon1 < 1) + length(lon);
            nbLon2 = indLon + nbDist + 1;
            nbLon2(nbLon2 > length(lon)) = nbLon2(nbLon2 > length(lon)) - length(lon);
            nbLon2(nbLon2 < 1) = nbLon2(nbLon2 < 1) + length(lon);
        
            Nedge = length(lonE);
            Uq1=zeros(Nedge,1);
            Uq2=zeros(Nedge,1);
            Uq3=zeros(Nedge,1);
            Uq4=zeros(Nedge,1);
            for i=1:Nedge
                Uq1(i)=UBOT(nbLon1(i),nbLat1(i));
                Uq2(i)=UBOT(nbLon2(i),nbLat1(i));
                Uq3(i)=UBOT(nbLon1(i),nbLat2(i));
                Uq4(i)=UBOT(nbLon2(i),nbLat2(i));
            end
            
            Vq1=zeros(Nedge,1);
            Vq2=zeros(Nedge,1);
            Vq3=zeros(Nedge,1);
            Vq4=zeros(Nedge,1);
            for i=1:Nedge
                Vq1(i)=VBOT(nbLon1(i),nbLat1(i));
                Vq2(i)=VBOT(nbLon2(i),nbLat1(i));
                Vq3(i)=VBOT(nbLon1(i),nbLat2(i));
                Vq4(i)=VBOT(nbLon2(i),nbLat2(i));
            end
            
            % map Vw to x, y and z in Euclidean
            el=lat(nbLat1);
            az=lon(nbLon1);
            cartVec1 = [cosd(el).*(-sind(az)).*Uq1-sind(el).*cosd(az).*Vq1, cosd(el).*cosd(az).*Uq1-sind(el).*sind(el).*Vq1, cosd(el).*Vq1]*pi/180;
            el=lat(nbLat1);
            az=lon(nbLon2);
            cartVec2 = [cosd(el).*(-sind(az)).*Uq2-sind(el).*cosd(az).*Vq2, cosd(el).*cosd(az).*Uq2-sind(el).*sind(el).*Vq2, cosd(el).*Vq2]*pi/180;
            el=lat(nbLat2);
            az=lon(nbLon1);
            cartVec3 = [cosd(el).*(-sind(az)).*Uq3-sind(el).*cosd(az).*Vq3, cosd(el).*cosd(az).*Uq3-sind(el).*sind(el).*Vq3, cosd(el).*Vq3]*pi/180;
            el=lat(nbLat2);
            az=lon(nbLon2);
            cartVec4 = [cosd(el).*(-sind(az)).*Uq4-sind(el).*cosd(az).*Vq4, cosd(el).*cosd(az).*Uq4-sind(el).*sind(el).*Vq4, cosd(el).*Vq4]*pi/180;

            % Set basis
            gridLandSeaE = repmat(obj.lMat, 1, 3)';
            gridLandSeaE = gridLandSeaE(:);
            obj.lMatE = gridLandSeaE;
            obj.QMaker.X1 = [cartVec(:,1).*(1-obj.lMatE), cartVec1(:,1).*(1-obj.lMatE), cartVec2(:,1).*(1-obj.lMatE), cartVec3(:,1).*(1-obj.lMatE), cartVec4(:,1).*(1-obj.lMatE)];
            obj.QMaker.X2 = [cartVec(:,2).*(1-obj.lMatE), cartVec1(:,2).*(1-obj.lMatE), cartVec2(:,2).*(1-obj.lMatE), cartVec3(:,2).*(1-obj.lMatE), cartVec4(:,2).*(1-obj.lMatE)];
            obj.QMaker.X3 = [cartVec(:,3).*(1-obj.lMatE), cartVec1(:,3).*(1-obj.lMatE), cartVec2(:,3).*(1-obj.lMatE), cartVec3(:,3).*(1-obj.lMatE), cartVec4(:,3).*(1-obj.lMatE)];
            
            obj.QMaker.nbNum = nbNum;
            obj.QMaker.nodeNum = nodeNum;
            % Set indicies
            currMax = max([obj.QMaker.idxRho, obj.QMaker.idxSigma, obj.QMaker.idxDVx, obj.QMaker.idxDVy, obj.QMaker.idxDVz]);
            obj.QMaker.idxNN = (currMax + (1 : (size(obj.QMaker.X1,2) * nodeNum + 2 * nodeNum + 1)*3));
        end
        
        function addWind(obj, sepVec)
           % ADDWIND
           %    Add wind covariate
           
            % Load zonal wind
            load('../Data/UBOT.mat');

            % Load meridional wind
            load('../Data/VBOT.mat');

            % Expand longitude and latitude to a mesh
            [gLat, gLon] = meshgrid(lat, lon);
            
            % Fix edge
            gLon = [gLon;gLon(1,:)+360];
            gLat = [gLat;gLat(1,:)];
            UBOT = [UBOT;UBOT(1,:)];
            VBOT = [VBOT;VBOT(1,:)];
            
            % Re-grid
            [az, el, ~] = cart2sph(obj.QMaker.meshELoc(:,1), obj.QMaker.meshELoc(:,2), obj.QMaker.meshELoc(:, 3));
            lonE = az*360/(2*pi);
            lonE(lonE < 0) = lonE(lonE < 0) + 360;
            latE = el*360/(2*pi);
            Uq = interp2(gLat, gLon, UBOT, latE, lonE);
            Vq = interp2(gLat, gLon, VBOT, latE, lonE);
            
            % Map to cartesian vector field
            cartVec = [cos(el).*(-sin(az)).*Uq-sin(el).*cos(az).*Vq, cos(el).*cos(az).*Uq-sin(el).*sin(el).*Vq, cos(el).*Vq]*pi/180;
           
            % Set basis
            gridLandSeaE = repmat(obj.lMat, 1, 3)';
            gridLandSeaE = gridLandSeaE(:);
            obj.lMatE = gridLandSeaE;
            if(~sepVec)
                obj.QMaker.basisDVx = [obj.QMaker.basisDVx, cartVec(:,1).*(1-obj.lMatE)];
                obj.QMaker.basisDVy = [obj.QMaker.basisDVy, cartVec(:,2).*(1-obj.lMatE)];
                obj.QMaker.basisDVz = [obj.QMaker.basisDVz, cartVec(:,3).*(1-obj.lMatE)];
            else
                zeroVec = 0*cartVec(:,1);
                obj.QMaker.basisDVx = [obj.QMaker.basisDVx, cartVec(:,1).*(1-obj.lMatE), zeroVec,                     zeroVec];
                obj.QMaker.basisDVy = [obj.QMaker.basisDVy, zeroVec,                     cartVec(:,2).*(1-obj.lMatE), zeroVec];
                obj.QMaker.basisDVz = [obj.QMaker.basisDVz, zeroVec,                     zeroVec,                     cartVec(:,3).*(1-obj.lMatE)];
            end
                
            
              
            % Set indicies
            currMax = max([obj.QMaker.idxRho, obj.QMaker.idxSigma, obj.QMaker.idxDVx, obj.QMaker.idxDVy, obj.QMaker.idxDVz]);
            if(~sepVec)
                obj.QMaker.idxDVx = [obj.QMaker.idxDVx, currMax+1];
                obj.QMaker.idxDVy = [obj.QMaker.idxDVy, currMax+1];
                obj.QMaker.idxDVz = [obj.QMaker.idxDVz, currMax+1];
            else
                obj.QMaker.idxDVx = [obj.QMaker.idxDVx, currMax+[1,2,3]];
                obj.QMaker.idxDVy = [obj.QMaker.idxDVy, currMax+[1,2,3]];
                obj.QMaker.idxDVz = [obj.QMaker.idxDVz, currMax+[1,2,3]];
            end
        end
        
        function addDistanceFromBorder(obj, sig)
           % ADDDISTANCEFROMBORDER
           %    Add distance from border as covariate
            %% Get distances from continent borders
            S = shaperead('../Data/continents/continent');
            X1 = S(1).X;
            Y1 = S(1).Y;
            for i = 2:7
                [xa, ya] = polybool('union', X1, Y1, S(i).X, S(i).Y);
                X1 = xa;
                Y1 = ya;
            end
            X1 = [X1,S(8).X];
            Y1 = [Y1,S(8).Y];
            
            % Remove line at latitude -90
            idx = (Y1<-90+1e-3);
            X1(idx) = [];
            Y1(idx) = [];
            

            cBorder = cell(1);
            for i = 1:1
                [x, y, z] = sph2cart(X1*pi/180, Y1*pi/180, 1);
                cBorder{i} = [x', y', z'];
            end
            
            % Get coordinates of centres of vertices
			xCenter = obj.QMaker.meshCLoc;

            % Find distances
            dGrid = abs(inf*xCenter(:,1));
            tic;
            for j = 1:length(dGrid)
                % Subsample
                idx = floor(linspace(1, size(cBorder{1}, 1), 10000));
                dGrid(j) = min(dGrid(j), min(pdist2([xCenter(j,1), xCenter(j,2), xCenter(j,3)], cBorder{1}(idx,:))));
            end
            [~, el, ~] = cart2sph(obj.QMaker.meshCLoc(:,1), obj.QMaker.meshCLoc(:,2), obj.QMaker.meshCLoc(:, 3));
            toc;

            bScale = exp(-dGrid.^2/sig^2);
         
            % Set basis
            obj.QMaker.basisRho = [obj.QMaker.basisRho, bScale];
            obj.QMaker.basisSigma = [obj.QMaker.basisSigma, bScale];
             
            % Set indicies
            currMax = max([obj.QMaker.idxRho, obj.QMaker.idxSigma, obj.QMaker.idxDVx, obj.QMaker.idxDVy, obj.QMaker.idxDVz]);
            obj.QMaker.idxRho = [obj.QMaker.idxRho, currMax+1];
            obj.QMaker.idxSigma = [obj.QMaker.idxSigma, currMax+2];
        end

		function calcBasis(obj, rhoMaxL, sigmaMaxL, vectorMaxL)
			% CALCBASIS
			%    Calculate basis functions for Rho and Sigma (currently)

			% Find az, el of centres
			[az, el, ~] = cart2sph(obj.QMaker.meshCLoc(:,1), obj.QMaker.meshCLoc(:,2), obj.QMaker.meshCLoc(:, 3));

			% Spherical harmonics!
			obj.QMaker.basisRho = [];
			for i = 0:rhoMaxL
				obj.QMaker.basisRho = [obj.QMaker.basisRho, SPDE.SphericalHarmonics.sphericalReal(i, az, pi/2-el)];
			end
			obj.QMaker.basisSigma = [];
			for i = 0:sigmaMaxL
				obj.QMaker.basisSigma = [obj.QMaker.basisSigma, SPDE.SphericalHarmonics.sphericalReal(i, az, pi/2-el)];
			end

			% Vector spherical harmonics!
			[az, el, ~] = cart2sph(obj.QMaker.meshELoc(:,1), obj.QMaker.meshELoc(:,2), obj.QMaker.meshELoc(:,3));
			obj.QMaker.basisDVx = [];
			obj.QMaker.basisDVy = [];
			obj.QMaker.basisDVz = [];
			if(vectorMaxL == 0)
				zer = zeros(size(obj.QMaker.basisRho(:,1),1)*3,1);
				obj.QMaker.basisDVx = [obj.QMaker.basisDVx, zer];
				obj.QMaker.basisDVy = [obj.QMaker.basisDVy, zer];
				obj.QMaker.basisDVz = [obj.QMaker.basisDVz, zer];
			else
				for i = 1:vectorMaxL
					[~, vx, vy, vz] = SPDE.SphericalHarmonics.sphericalR3(i, az, pi/2-el);
					obj.QMaker.basisDVx = [obj.QMaker.basisDVx, vx];
					obj.QMaker.basisDVy = [obj.QMaker.basisDVy, vy];
					obj.QMaker.basisDVz = [obj.QMaker.basisDVz, vz];
				end
			end

			% Indicies!
			obj.QMaker.idxRho = 1:size(obj.QMaker.basisRho, 2);
			obj.QMaker.idxSigma = obj.QMaker.idxRho(end)+(1:size(obj.QMaker.basisSigma, 2));
			obj.QMaker.idxDVx = obj.QMaker.idxSigma(end)+(1:size(obj.QMaker.basisDVx, 2));
			obj.QMaker.idxDVy = obj.QMaker.idxDVx;
			obj.QMaker.idxDVz = obj.QMaker.idxDVx;
        end
        
        function [A] = getObsMatrix(obj, pos)
			% GETOBSMATRIX
            %   get observation matrix for mapping to pos

			% Make map from vertex to triangles
			nT = size(obj.QMaker.meshTV, 1);
			nV = size(obj.QMaker.meshVLoc, 1);
			vt = [obj.QMaker.meshTV(:,1), (1:nT)'; obj.QMaker.meshTV(:,2), (1:nT)'; obj.QMaker.meshTV(:,3), (1:nT)'];
			vtTmp = unique(vt, 'rows');
			vt = cell(nV, 1);
			idxPrev = 1;
			idxNext = 1;
			vtTmp = [vtTmp; -1, -1];
			for i = 1:nV
				while(vtTmp(idxNext,1) == vtTmp(idxPrev,1))
					idxNext = idxNext+1;
				end
				vt{i} = vtTmp(idxPrev:(idxNext-1), 2);
				idxPrev = idxNext;
			end

			% Find closest vertex for each observation position
			idx = knnsearch(obj.QMaker.meshVLoc, pos, 'K', 50);

			% For each observation check the triangles connected to the closest vertex
			trClose = zeros(length(idx), 1);
			tic;
			for i = 1:size(idx,1)
				minDist = inf;
				minTr = -1;
                for j = 1:size(idx,2)
                    for tr = vt{idx(i,j)}'
                        p1 = obj.QMaker.meshVLoc(obj.QMaker.meshTV(tr,1),:);
                        p2 = obj.QMaker.meshVLoc(obj.QMaker.meshTV(tr,2),:);
                        p3 = obj.QMaker.meshVLoc(obj.QMaker.meshTV(tr,3),:);
                        d = SPDE.pointTriangleDistance([p1; p2; p3], pos(i,:));
                        if(d < minDist)
                            minDist = d;
                            minTr = tr;
                        end
                    end
                    trClose(i) = minTr;
                end
			end
			toc;
            
            % Assemble matrix
            nObs = size(pos, 1);
            A = sparse(1:nObs, trClose, 1, nObs, nT);
        end
        
        function [Slin] = setObservation(obj, pos, obs)
            % Make map from vertex to triangles
			nT = size(obj.QMaker.meshTV, 1);
			nV = size(obj.QMaker.meshVLoc, 1);
			vt = [obj.QMaker.meshTV(:,1), (1:nT)'; obj.QMaker.meshTV(:,2), (1:nT)'; obj.QMaker.meshTV(:,3), (1:nT)'];
			vtTmp = unique(vt, 'rows');
			vt = cell(nV, 1);
			idxPrev = 1;
			idxNext = 1;
			vtTmp = [vtTmp; -1, -1];
			for i = 1:nV
				while(vtTmp(idxNext,1) == vtTmp(idxPrev,1))
					idxNext = idxNext+1;
				end
				vt{i} = vtTmp(idxPrev:(idxNext-1), 2);
				idxPrev = idxNext;
			end

			% Find closest vertex for each observation position
			idx = knnsearch(obj.QMaker.meshVLoc, pos, 'K', 50);

			% For each observation check the triangles connected to the closest vertex
			trClose = zeros(length(idx), 1);
			tic;
			for i = 1:size(idx,1)
				minDist = inf;
				minTr = -1;
                for j = 1:size(idx,2)
                    for tr = vt{idx(i,j)}'
                        p1 = obj.QMaker.meshVLoc(obj.QMaker.meshTV(tr,1),:);
                        p2 = obj.QMaker.meshVLoc(obj.QMaker.meshTV(tr,2),:);
                        p3 = obj.QMaker.meshVLoc(obj.QMaker.meshTV(tr,3),:);
                        d = SPDE.pointTriangleDistance([p1; p2; p3], pos(i,:));
                        if(d < minDist)
                            minDist = d;
                            minTr = tr;
                        end
                    end
                    trClose(i) = minTr;
                end
			end
			toc;
            
            % Calculate piece-wise linear S
            nLoc = length(trClose);
            tmpI = zeros(nLoc*3,1);
            tmpJ = tmpI;
            tmpV = tmpJ;
            vLoc = obj.QMaker.meshVLoc;
            tv = obj.QMaker.meshTV;

            for i = 1:nLoc
                % Which cell
                cIdx = trClose(i);

                % Get vertexes
                v1 = vLoc(tv(cIdx,1),:);
                v2 = vLoc(tv(cIdx,2),:);
                v3 = vLoc(tv(cIdx,3),:);

                % Shift coordinate system
                vec1 = (v2-v1);
                vec2 = (v3-v1);

                % Calculate barycentric coordinates
                cord = abs([vec1',vec2']\(pos(i,:)'-v1'));
                cord = [1-sum(cord);cord];

                % create linear combination needed
                ttJ = [];
                ttV = [];
                for j = 1:3
                    ttJ = [ttJ; tv(cIdx,j)];
                    ttV = [ttV; cord(j)];
                end
                ttI = ones(size(ttJ,1),1)*i;

                % Update full vectors
                tmpI(3*(i-1)+(1:3)) = ttI;
                tmpJ(3*(i-1)+(1:3)) = ttJ;
                tmpV(3*(i-1)+(1:3)) = ttV;
            end
            
            % Assemble matrix
            Slin = sparse(tmpI,tmpJ,tmpV, nLoc, nV);
            
            % Store each observation
            obj.A = Slin*obj.QMaker.AtoV;
            obj.obs = cell(size(obs,2),1);
            for i = 1:size(obs, 2)
                % Store observations
                obj.obs{i} = obs(:,i);
            end
            
            % TODO: temporaraily revert to old way of dealing with things
            obj.A = sparse(1:nLoc, trClose, ones(nLoc, 1), nLoc, nT);
        end
        
        function [Q, S, Qc, Qn] = combMatrix(obj, par, X, tauMu, shouldCor)
            % Remove parameters in variance
                par(obj.QMaker.idxSigma) = 0;
                
                % Get current precision matrices
                Q = obj.QMaker.makeQ(par(1:end-2));
                     
                % Scale
                Qinv = obj.safeQinv(Q);
                Dc = spdiags(sqrt(diag(Qinv)), 0, size(Qinv,1), size(Qinv,2));
                
                % Correct precision matrices
                Q = Dc*Q*Dc;
                
                % Get noise precision
                tauNsea = exp(par(end-1));
                tauNland = exp(par(end));
                numObs = size(obj.A, 1);
                [i,j,~] = find(obj.A);
                [~,I] = sort(i);
                lPix = obj.lMat(j(I));
                sPix = 1-lPix;
                QnSea = spdiags(sPix, 0, numObs, numObs)*tauNsea;
                QnLand = spdiags(lPix, 0, numObs, numObs)*tauNland;
                %diag(QnSea) = obj.A(:,tSea)*tauNsea;
                %diag(QnLand) = tLand*tauNland;
                Qn = QnSea + QnLand;
                
                % Set prior precision of latent variables
                Qmu = tauMu*speye(size(X,2));

                % Include covariates in latent model
                nQ = size(Q, 1);
                nCov = size(Qmu, 1);
                Q = [Q, sparse(nQ, nCov);
                    sparse(nCov, nQ), Qmu];
                
                % Total selection matrix
                S = [obj.A, X];

                % Update with information about observations
                Qc = Q + S'*Qn*S;
        end

		function [val, gradLL] = logLikelihood(obj, par, X, tauMu, rmIdx, shouldCor, h, useCentral)
            diary off
            diary on
			% LOGLIKELIHOOD
			%    Calculate minus log-likelihood based on parameters
            %% Shared for all realizations
                
                % Remove parameters in variance
                par(obj.QMaker.idxSigma) = 0;
                
                % Get current precision matrices
                Q = obj.QMaker.makeQ(par(1:end-2));
                     
                % Scale
                try
                    Qinv = SPDE.INLA.sparseQinv(Q);
                catch
                    warning('Q not invertible logLikelihood');
                    val = inf;
                    gradLL = eye(length(par), 1);
                    return;
                end
                Dc = spdiags(sqrt(diag(Qinv)), 0, size(Qinv,1), size(Qinv,2));
                
                % Correct precision matrices
                Q = Dc*Q*Dc;
                
                % Get noise precision
                tauNsea = exp(par(end-1));
                tauNland = exp(par(end));
                numObs = size(obj.A, 1);
                [i,j,~] = find(obj.A);
                [~,I] = sort(i);
                lPix = obj.lMat(j(I));
                sPix = 1-lPix;
                QnSea = spdiags(sPix, 0, numObs, numObs)*tauNsea;
                QnLand = spdiags(lPix, 0, numObs, numObs)*tauNland;
                %diag(QnSea) = obj.A(:,tSea)*tauNsea;
                %diag(QnLand) = tLand*tauNland;
                Qn = QnSea + QnLand;
                
                % Set prior precision of latent variables
                Qmu = tauMu*speye(size(X,2));

                % Include covariates in latent model
                nQ = size(Q, 1);
                nCov = size(Qmu, 1);
                Q = [Q, sparse(nQ, nCov);
                    sparse(nCov, nQ), Qmu];
                
                % Total selection matrix
                S = [obj.A, X];

                % Update with information about observations
                Qc = Q + S'*Qn*S;
                
                % Calculate Cholesky factorizations
                [L, Lc, p] = obj.safeChol(Q, Qc);
                Ln = chol(Qn);

           %% Separate for each realization
                % Iterate through each realization
                Lct = Lc';
                St = S';
                val = length(obj.obs)*(sum(log(diag(L)))+sum(log(diag(Ln)))-sum(log(diag(Lc))));
                valTmp = ones(length(obj.obs),1);
                tmpObs = obj.obs;
                parfor i = 1:length(obj.obs)
                    % Calculate conditional mean
                    muIy = Lct\(Lc\(St(p,:)*Qn*tmpObs{i}));
                    muIy(p) = muIy;
                    muAll{i} = muIy;

                    % Calculate exponent term
                    valTmp(i) = - 0.5*muIy'*Q*muIy - 0.5*(tmpObs{i}-S*muIy)'*Qn*(tmpObs{i}-S*muIy);
                end
                val = val + sum(valTmp);

                
            %% Calculate negative log-likelihood
			val = -val;
            
            % Calculate per observation
            val = val/length(obj.obs);
            
			%% Find gradient if needed
			if(nargout > 1)        
                % Precompute inverses
                [Qinv, Qcinv] = obj.safeQinv(Q, Qc);
                
				% Initialize storage
                nPar = length(par);
				gradLL = zeros(length(par), 1);
                %tic;
                %% Sea
                    % Derivatives w.r.t. \rho
                    tmpGrad = zeros(length(obj.QMaker.idxRho), 1);
                    parfor i = 1:length(obj.QMaker.idxRho)
                        % Get derivative of Q ignoring scale
                        QderTmp = obj.QMaker.rhoDer(par(obj.QMaker.idxRho), i);
                        Qder = Dc*QderTmp*Dc;

                        % Derivative of scaling
                        parSea2 = par;
                        parSea2(obj.QMaker.idxSigma) = 0;
                        parSea2(obj.QMaker.idxRho(i)) = parSea2(obj.QMaker.idxRho(i)) + h;
                        Qsea2 = obj.QMaker.makeQ(parSea2);
                        Qinv2 = obj.safeQinv(Qsea2);
                        Dc1H = spdiags(sqrt(diag(Qinv2)), 0, size(Qinv2,1), size(Qinv2,2));
                        deltaD = (Dc1H-Dc)/h;
                        if(useCentral)
                            parSea2(obj.QMaker.idxRho(i)) = parSea2(obj.QMaker.idxRho(i)) - 2*h;
                            Qsea2 = obj.QMaker.makeQ(parSea2);
                            Qinv2 = obj.safeQinv(Qsea2);
                            Dcm1H = spdiags(sqrt(diag(Qinv2)), 0, size(Qinv2,1), size(Qinv2,2));

                            deltaD = (Dc1H-Dcm1H)/(2*h);
                        end

                        Q1 = obj.QMaker.makeQ(par);
                        Qder2 = deltaD*Q1*Dc;
                        Qder2 = Qder2 + Qder2';

                        % Full derivative 
                        Qder = Qder + Qder2;
                        Qder = (Qder+Qder')/2;

%                         % Do derivative given Qder for each realization
%                         for j = 1:length(obj.obs)
%                             tmpGrad(i) = tmpGrad(i) + obj.diffWRTQpar(Qder, Qinv, Qcinv, muAll{j});
%                         end
%                         
                        
                        % Derivative of determinants
                        tmpValGrad = length(obj.obs)*0.5*sum(sum((Qinv-Qcinv).*Qder));
                        for j = 1:length(obj.obs)
                            tmpValGrad = tmpValGrad - 0.5*muAll{j}'*Qder*muAll{j};
                        end
                        tmpGrad(i) = tmpValGrad;
                    end
                    gradLL(obj.QMaker.idxRho) = tmpGrad;
                    
                    % Derivative w.r.t. ani
                    tmpGrad = zeros(length(obj.QMaker.idxDVx),1);
                    parfor i = 1:length(obj.QMaker.idxDVx)
                        % Get derivative of Q
                        QderTmp = obj.QMaker.aniDer(par(obj.QMaker.idxDVx), i);
                        Qder = Dc*QderTmp*Dc;

                        % Derivative of scaling
                        parSea2 = par;
                        parSea2(obj.QMaker.idxSigma) = 0;
                        parSea2(obj.QMaker.idxDVx(i)) = parSea2(obj.QMaker.idxDVx(i)) + h;
                        Qsea2 = obj.QMaker.makeQ(parSea2);
                        Qinv2 = obj.safeQinv(Qsea2);
                        Dc1H = spdiags(sqrt(diag(Qinv2)), 0, size(Qinv2,1), size(Qinv2,2));
                        deltaD = (Dc1H-Dc)/h;
                        if(useCentral)
                            parSea2(obj.QMaker.idxDVx(i)) = parSea2(obj.QMaker.idxDVx(i)) - 2*h;
                            Qsea2 = obj.QMaker.makeQ(parSea2);
                            Qinv2 = obj.safeQinv(Qsea2);
                            Dcm1H = spdiags(sqrt(diag(Qinv2)), 0, size(Qinv2,1), size(Qinv2,2));

                            deltaD = (Dc1H-Dcm1H)/(2*h);
                        end

                        Q1 = obj.QMaker.makeQ(par);
                        Qder2 = deltaD*Q1*Dc;
                        Qder2 = Qder2 + Qder2';

                        % Full derivative 
                        Qder = Qder + Qder2;
                        Qder = (Qder+Qder')/2;
                        
%                         % Do derivative given Qder for each realization
%                         for j = 1:length(obj.obs)
%                             tmpGrad(i) = tmpGrad(i) + obj.diffWRTQpar(Qder, Qinv, Qcinv, muAll{j});
%                         end
%                         
                        % Derivative of determinants
                        tmpValGrad = length(obj.obs)*0.5*sum(sum((Qinv-Qcinv).*Qder));
                        for j = 1:length(obj.obs)
                            tmpValGrad = tmpValGrad - 0.5*muAll{j}'*Qder*muAll{j};
                        end
                        tmpGrad(i) = tmpValGrad;
                    end
                    gradLL(obj.QMaker.idxDVx) = tmpGrad;

                    %% Noise    
                    % Derivative w.r.t. noise
                    QnDer = QnSea;
                    QnoiseInv = speye(numObs)/tauNsea;
%                     tic;
%                     for j = 1:length(obj.obs)
%                         gradLL(end-1) = gradLL(end-1) + obj.diffWRTNoisePar(QnDer, QnoiseInv, Qcinv, muAll{j}, S, obj.obs{j});
%                     end
%                     toc;
%                     
                    tmpValGrad = length(obj.obs)*(0.5*sum(sum(QnoiseInv.*QnDer)) - 0.5*sum(sum(Qcinv.*(S'*QnDer*S))));
                    for j = 1:length(obj.obs)
                        tmpValGrad = tmpValGrad - 0.5*(obj.obs{j}-S*muAll{j})'*QnDer*(obj.obs{j}-S*muAll{j});
                    end
                    gradLL(end-1) = tmpValGrad;

                    QnDer = QnLand;
                    QnoiseInv = speye(numObs)/tauNland;
%                     for j = 1:length(obj.obs)
%                         gradLL(end) = gradLL(end) + obj.diffWRTNoisePar(QnDer, QnoiseInv, Qcinv, muAll{j}, S, obj.obs{j});
%                     end
                    tmpValGrad = length(obj.obs)*(0.5*sum(sum(QnoiseInv.*QnDer)) - 0.5*sum(sum(Qcinv.*(S'*QnDer*S))));
                    for j = 1:length(obj.obs)
                        tmpValGrad = tmpValGrad - 0.5*(obj.obs{j}-S*muAll{j})'*QnDer*(obj.obs{j}-S*muAll{j});
                    end
                    gradLL(end) = tmpValGrad;
                   
                   %disp('gradients:')
                   %toc;
                %% The rest
                    % Change sign
                    gradLL = -gradLL;
                    
                    % Calculate per observation
                    gradLL = gradLL/length(obj.obs);
                    
                    %display('Gradient')
                    % gradLL
            
                %% Delete fixed parameters
                gradLL(rmIdx) = [];
            end	
        end
        
        function [pen, penGrad] = penalty(obj, par, lam)
            % PENALTY
            %    Add penalty on standard deviation. We don't
            %    want intrinsic models
            
            % Calculate approx sigma
            theta = par(obj.QMaker.idxSigma);
            DSig = obj.QMaker.makeDSigma(theta, 0);
            sig = diag(DSig);
            
            % Get weights
            DV = obj.QMaker.makeDV(0);
            w = diag(DV);
            
            % average standard deviation
            sigAv = sum(w.*sig)/sum(w);
            
            % Calculate penalty
            pen = lam*sigAv;
            
            % Calculate gradient if necessary
            grad = zeros(length(theta), 1);
   			if(nargout > 1)
                for i = 1:length(theta)
                    % Derivative of sigma 
                    sigDer = obj.QMaker.basisRho(:,i).*sig;
                    sigAvDer = sum(w.*sigDer)/sum(w);

                    % Gradient
                    grad(i) = lam*sigAvDer;
                end
            end
            
            penGrad = grad;
        end

		function [val] = diffWRTQpar(obj, Qder, Qinv, Qcinv, muIy)
			% DIFFWRTQPAR
			%    Derivative of log-likelihood given derivative of Q

			% Derivative of determinants
			nonZer = size(Qder,1);
			val = 0.5*sum(sum((Qinv(1:nonZer,1:nonZer)-Qcinv(1:nonZer,1:nonZer)).*Qder));

			% Derivative of quadratic forms
			val = val - 0.5*muIy(1:nonZer)'*Qder*muIy(1:nonZer);
		end

		function [val] = diffWRTNoisePar(obj, QnDer, QnoiseInv, Qcinv, muIy, S, y)
			% Derivative of determinants
			val = 0.5*sum(sum(QnoiseInv.*QnDer)) - 0.5*sum(sum(Qcinv.*(S'*QnDer*S)));

			% Derivative of quadratic forms
			val = val - 0.5*(y-S*muIy)'*QnDer*(y-S*muIy);
        end
	end

	methods(Static)
        function OptNStat = makeNonStatModel(vLoc, tt, tv, loc, allObs, barrier, tensorPar, orderSpec)
            % Need this for simulation
            OptNStat = SPDE.Optimizer(vLoc, tt, tv, tensorPar);

            % Calculate basis
            if(nargin < 8)
                OptNStat.calcBasis(4,0,4);
            else
                OptNStat.calcBasis(orderSpec(1), orderSpec(2), orderSpec(3));
            end
            
            % Add land/sea
            tmpBas1 = OptNStat.QMaker.basisRho;
            tmpBas2 = tmpBas1;
            if(~barrier)
                tmpBas1 = tmpBas1(:,1);
            end
            for i = 1:size(tmpBas1,2)
                tmpBas1(:,i) = tmpBas1(:,i).*OptNStat.lMat;
            end
            for i = 1:size(tmpBas2,2)
                tmpBas2(:,i) = tmpBas2(:,i).*(1-OptNStat.lMat);
            end
            OptNStat.QMaker.basisRho = [tmpBas1, tmpBas2];
            cMax = max(OptNStat.QMaker.idxDVx)+1;
            OptNStat.QMaker.idxRho = [OptNStat.QMaker.idxRho, cMax:(cMax+size(tmpBas1,2)-1)];
            
            
            tmpBas1 = OptNStat.QMaker.basisDVx;
            tmpBas2 = tmpBas1;
            if(~barrier)
                tmpBas1 = zeros(size(tmpBas1, 1), 1);
            end
            for i = 1:size(tmpBas1,2)
                tmpBas1(:,i) = tmpBas1(:,i).*OptNStat.lMatE;
            end
            for i = 1:size(tmpBas2,2)
                tmpBas2(:,i) = tmpBas2(:,i).*(1-OptNStat.lMatE);
            end
            OptNStat.QMaker.basisDVx = [tmpBas1, tmpBas2];

            tmpBas1 = OptNStat.QMaker.basisDVy;
            tmpBas2 = tmpBas1;
            if(~barrier)
                tmpBas1 = zeros(size(tmpBas1, 1), 1);
            end
            for i = 1:size(tmpBas1,2)
                tmpBas1(:,i) = tmpBas1(:,i).*OptNStat.lMatE;
            end
            for i = 1:size(tmpBas2,2)
                tmpBas2(:,i) = tmpBas2(:,i).*(1-OptNStat.lMatE);
            end 
            OptNStat.QMaker.basisDVy = [tmpBas1, tmpBas2];
            
            tmpBas1 = OptNStat.QMaker.basisDVz;
            tmpBas2 = tmpBas1;
            if(~barrier)
                tmpBas1 = zeros(size(tmpBas1, 1), 1);
            end
            for i = 1:size(tmpBas1,2)
                tmpBas1(:,i) = tmpBas1(:,i).*OptNStat.lMatE;
            end
            for i = 1:size(tmpBas2,2)
                tmpBas2(:,i) = tmpBas2(:,i).*(1-OptNStat.lMatE);
            end
            OptNStat.QMaker.basisDVz = [tmpBas1, tmpBas2];
            
            cMax = max(OptNStat.QMaker.idxRho)+1;
            OptNStat.QMaker.idxDVx = [OptNStat.QMaker.idxDVx, cMax:(cMax+size(tmpBas1,2)-1)];
            OptNStat.QMaker.idxDVy = OptNStat.QMaker.idxDVx;
            OptNStat.QMaker.idxDVz = OptNStat.QMaker.idxDVx;

            % Add barrier
            if(barrier)
                bSize = 0.03;
                bIdx = OptNStat.dMat <= bSize;
                OptNStat.QMaker.basisRho = [OptNStat.QMaker.basisRho, bIdx];
                OptNStat.QMaker.idxRho = [OptNStat.QMaker.idxRho, 148];
            end

            % Set observations (SLOW)
            OptNStat.setObservation(loc, allObs);

        end
        
        function OptStat = makeStationaryModel(vLoc, tt, tv, loc, allObs, barrier, tensorPar)
            % Use general constructor
            OptStat = SPDE.Optimizer(vLoc, tt, tv, tensorPar);

            % Calculate basis
            OptStat.calcBasis(0,0,0);

            % Set observations (SLOW)
            OptStat.setObservation(loc, allObs);

            % Add land/sea covariate
            OptStat.addLandSeaCov();

            % Add barrier
            if(barrier)
                bSize = 0.03;
                bIdx = OptStat.dMat <= bSize;
                OptStat.QMaker.basisRho = [OptStat.QMaker.basisRho, bIdx];
                OptStat.QMaker.idxRho = [OptStat.QMaker.idxRho, 5];
            end
        end

        function OptNN = makeNN(vLoc, tt, tv, loc, allObs, barrier, tensorPar, orderSpec)
            % Use general constructor
            OptNN = SPDE.Optimizer(vLoc, tt, tv, tensorPar);

            % Just use intercept for each parameter
            OptNN.QMaker = SPDE.PrecisionMaker_NN2(vLoc, tt, tv, tensorPar);
            
            % Calculate basis
            if(nargin < 8)
                OptNN.calcBasis(4,0,4);
            else
                OptNN.calcBasis(orderSpec(1), orderSpec(2), orderSpec(3));
            end

            % Add land/sea
            tmpBas1 = OptNN.QMaker.basisRho;
            tmpBas2 = tmpBas1;
            if(~barrier)
                tmpBas1 = tmpBas1(:,1);
            end
            for i = 1:size(tmpBas1,2)
                tmpBas1(:,i) = tmpBas1(:,i).*OptNN.lMat;
            end
            for i = 1:size(tmpBas2,2)
                tmpBas2(:,i) = tmpBas2(:,i).*(1-OptNN.lMat);
            end
            OptNN.QMaker.basisRho = [tmpBas1, tmpBas2];
            cMax = max(OptNN.QMaker.idxDVx)+1;
            OptNN.QMaker.idxRho = [OptNN.QMaker.idxRho, cMax:(cMax+size(tmpBas1,2)-1)];
            
            
            tmpBas1 = OptNN.QMaker.basisDVx;
            tmpBas2 = tmpBas1;
            if(~barrier)
                tmpBas1 = zeros(size(tmpBas1, 1), 1);
            end
            for i = 1:size(tmpBas1,2)
                tmpBas1(:,i) = tmpBas1(:,i).*OptNN.lMatE;
            end
            for i = 1:size(tmpBas2,2)
                tmpBas2(:,i) = tmpBas2(:,i).*(1-OptNN.lMatE);
            end
            OptNN.QMaker.basisDVx = [tmpBas1, tmpBas2];

            tmpBas1 = OptNN.QMaker.basisDVy;
            tmpBas2 = tmpBas1;
            if(~barrier)
                tmpBas1 = zeros(size(tmpBas1, 1), 1);
            end
            for i = 1:size(tmpBas1,2)
                tmpBas1(:,i) = tmpBas1(:,i).*OptNN.lMatE;
            end
            for i = 1:size(tmpBas2,2)
                tmpBas2(:,i) = tmpBas2(:,i).*(1-OptNN.lMatE);
            end 
            OptNN.QMaker.basisDVy = [tmpBas1, tmpBas2];
            
            tmpBas1 = OptNN.QMaker.basisDVz;
            tmpBas2 = tmpBas1;
            if(~barrier)
                tmpBas1 = zeros(size(tmpBas1, 1), 1);
            end
            for i = 1:size(tmpBas1,2)
                tmpBas1(:,i) = tmpBas1(:,i).*OptNN.lMatE;
            end
            for i = 1:size(tmpBas2,2)
                tmpBas2(:,i) = tmpBas2(:,i).*(1-OptNN.lMatE);
            end
            OptNN.QMaker.basisDVz = [tmpBas1, tmpBas2];
            
            cMax = max(OptNN.QMaker.idxRho)+1;
            OptNN.QMaker.idxDVx = [OptNN.QMaker.idxDVx, cMax:(cMax+size(tmpBas1,2)-1)];
            OptNN.QMaker.idxDVy = OptNN.QMaker.idxDVx;
            OptNN.QMaker.idxDVz = OptNN.QMaker.idxDVx;
            
            % Set observations (SLOW)
            OptNN.setObservation(loc, allObs);
        end
        
		function [isIn] = isInsideTriangle(p1, p2, p3, p0)
			% DPOINTTRIANGLE
			%    Calculate distance from triangle defined by p1, p2 and p3 to p0

			% Local axes
			e1 = (p2-p1);
			e2 = (p3-p1);

			% Local coordinates
			x = dot(e1, p0-p1);
			y = dot(e2, p0-p1);

			% Check if "inside"
			if( (x >= 0) & (x <= 1) & (y >= 0 & y <= 1))
				isIn = true;
			else
				isIn = false;
            end 
		end
			
		function [Qinv, Qcinv] = safeQinv(Q, Qc)
			try
				Qinv = SPDE.INLA.sparseQinv(Q);
			catch err
				display('Terrible things happened with Qinv = sparseQinv, trying again');
				try
					Qinv = SPDE.INLA.sparseQinv(Q);
				catch err2
					display('Failed again, just give som crazy value');
					Qinv = 0*Q;
				end
            end
            
            if(nargin > 1)
                try
                    Qcinv = SPDE.INLA.sparseQinv(Qc);
                catch err
                    display('Terrible things happened with Qcinv = sparseQinv, trying again');
                    try
                        Qcinv = SPDE.INLA.sparseQinv(Qc);
                    catch err2
                        display('Failed again, just give som crazy value');
                        Qcinv = 0*Qc;
                    end
                end
            end
        end
        
        function [Qinv, Qcinv] = safeQinvPar(Q, Qc)
            Qs = cell(2,1);
            Qs{1} = Q;
            Qs{2} = Qc;
            QQ = cell(2,1);
            parfor i = 1:2
                Qinv = Qs{i};
                try
                    Qinv = SPDE.INLA.sparseQinv(Qs{i});
                catch err
                    display('Terrible things happened with Qinv = sparseQinv, trying again');
                    try
                        Qinv = SPDE.INLA.sparseQinv(Qs{i});
                    catch err2
                        display('Failed again, just give som crazy value');
                        Qinv = 0*Qs{i};
                    end
                end
                QQ{i} = Qinv;
            end

			Qinv = QQ{1};
            Qcinv = QQ{2};
		end

		function [L, Lc, p] = safeChol(Q, Qc)
			% Permutation
			p = amd(Q);

			% Do stuff
			try
				L = chol(Q(p,p),'lower');
			catch err
				display('Could not do Cholesky factorization, trying again with extra \eps added')
				L = chol(Q(p,p)+max(max(Q))*1e-10*speye(size(Q)),'lower');
			end

			% Find necessary parts of Qinv, (Is this already calculated in call to INLA???)
			try
				Lc = chol(Qc(p,p), 'lower');
			catch err
				display('Could not do Cholesky factorization, trying again with extra \eps added');
				Lc = chol(Qc(p,p)+max(max(Qc))*1e-10*speye(size(Qc)),'lower');
			end
        end
        
        function [L, Lc, p] = safeCholPar(Q, Qc)
			% Permutation
			p = amd(Q);

			% Do stuff
            Qs = cell(2,1);
            Qs{1} = Q;
            Qs{2} = Qc;
            Ls = cell(2,1);
            parfor i = 1:2
                L = Qs{i};
                try
                    L = chol(Qs{i}(p,p),'lower');
                catch err
                    display('Could not do Cholesky factorization, trying again with extra \eps added')
                    L = chol(Qs{i}(p,p)+max(max(Q))*1e-10*speye(size(Q)),'lower');
                end
                Ls{i} = L;
            end
            L = Ls{1};
            Lc = Ls{2};
		end
	end
end
