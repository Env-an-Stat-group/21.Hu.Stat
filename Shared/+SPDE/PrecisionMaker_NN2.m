classdef PrecisionMaker_NN2 < handle
	% PRECISIONMAKER_NN
	%    New and fancy NN parametrization of non-stationarity

	properties
        % Type of parametrization for diffusion tensor
        tensorPar;
        
		% Geometry
		meshTV;
		meshVLoc;
		meshTT;
        meshVT;
		meshArea;
		meshCLoc;
		meshELoc;
        
        % Bases
		basisRho;
		basisSigma;
		basisDVx;
		basisDVy;
		basisDVz;

		% Indicies
		idxRho;
		idxSigma;
		idxDVx;
		idxDVy;
		idxDVz;

		% Store stuff to speed up computations
		DV;
		DVInv;
		DRho2;
		DRho2Inv;
		AH;
		DSigmaInv;
		Aop;
        
        % Neural net
        actFun;
        X1;
        X2;
        X3;
        idxNN;
        nbNum;
        nodeNum;
        
        % Vertex reconstruction matrices
        rIdx;
        rVal;
        startIdx;
        AtoV;
	end

	methods
		function [obj] = PrecisionMaker_NN2(vLoc, tt, tv, tensorPar)
			% PRECISIONMAKER
			%    Constructor based on Delenauy triangulation

            % Parametrization
            obj.tensorPar = tensorPar;
            
			% Store triangulation
			obj.meshTV = tv;
			obj.meshVLoc = vLoc;
			obj.meshTT = tt;

			% Calculate areas
			obj.meshArea = zeros(size(tv,1),1);
			for i = 1:length(obj.meshArea)
				p1 = vLoc(tv(i,1),:);
				p2 = vLoc(tv(i,2),:);
				p3 = vLoc(tv(i,3),:);
				obj.meshArea(i) = obj.triangleArea(p1, p2, p3);
			end

			% Calculate locations of centroids
			obj.meshCLoc = (vLoc(tv(:,1),:)+vLoc(tv(:,2),:)+vLoc(tv(:,3),:))/3.0;

			% Calculate locations of centroids of each edge
			eLoc1 = (vLoc(tv(:,2),:)+vLoc(tv(:,3),:))/2;
			eLoc2 = (vLoc(tv(:,1),:)+vLoc(tv(:,3),:))/2;
			eLoc3 = (vLoc(tv(:,1),:)+vLoc(tv(:,2),:))/2;
			obj.meshELoc = zeros(size(eLoc1,1)*3,size(eLoc1,2));
			obj.meshELoc(3*(1:size(eLoc1,1))-2,:) = eLoc1;
			obj.meshELoc(3*(1:size(eLoc1,1))-1,:) = eLoc2;
			obj.meshELoc(3*(1:size(eLoc1,1)),:)   = eLoc3;
            
            % Calculate vertex to triangles map
            obj.meshVT =  cell(max(max(tv)),1);
            for vt = 1:length(obj.meshVT)
                obj.meshVT{vt} = [];
            end
            for tr = 1:size(obj.meshTV,1)
                for j = 1:3
                    obj.meshVT{obj.meshTV(tr,j)} = [obj.meshVT{obj.meshTV(tr,j)};tr];
                end
            end
            
            % Calculate vertex reconstruction map
            rIdx = [];
            rVal = [];
            vIdx = [];
            startIdx = [1];
            for vt = 1:size(obj.meshVT,1)
                % Triangles that includes vertex
                tIdx = obj.meshVT{vt};
                
                % Weighting
                lam = obj.meshArea(tIdx);
                lam = lam/sum(lam);
                Lam = sparse(1:length(lam),1:length(lam),lam, length(lam), length(lam));
                
                % New coordinate system
                nVec = vLoc(vt,:)';
                v1 = zeros(3,1);
                v2 = zeros(3,1);
                if(abs(nVec(3)) ~= 1)
                    v1([1;2]) = [nVec(2)*-1, nVec(1)];
                else
                    v1(1) = 1;
                end
                v1 = v1/norm(v1);
                v2 = cross(nVec, v1);
                v2 = v2/norm(v2);
                
                % Project coordinates
                vOld = obj.meshCLoc(tIdx,:);
                for i = 1:length(tIdx)
                    vOld(i,:) = vOld(i,:)-nVec';
                end
                xNew = vOld*v1;
                yNew = vOld*v2;
                
                % Solve weighted least squares
                A = [ones(length(tIdx),1), xNew, yNew];
                weights = (A'*Lam*A)\(A'*Lam);
                
                rVal = [rVal; weights(1,:)'];
                rIdx = [rIdx; tIdx];
                vIdx = [vIdx; vt*ones(length(tIdx),1)];
                startIdx = [startIdx; startIdx(end)+length(tIdx)];
            end
            
            obj.actFun=@tanh;
            
            
            % Store values
            obj.rVal = rVal;
            obj.rIdx = rIdx;
            obj.AtoV = sparse(vIdx, rIdx, rVal, max(vIdx), max(rIdx));
            
            % Shift to zero-indexing
            obj.startIdx = startIdx-1;
            obj.rIdx = obj.rIdx-1;
        end

        function [vx, vy, vz] = calcCoeff(obj, parHM, parNN)
           
           p = obj.nbNum + 1;
           d = obj.nodeNum;
           n = size(obj.X1,1);
           parNum = p*d + d + d + 1;% + 1;
           
           % Calculate vx
           xWHidden = parNN(1 : p * d);
           xWHidden = reshape(xWHidden, p, d);
           xbHidden = parNN(p*d + (1 : d));
           xbHidden = repmat(xbHidden',n,1);
           xWOut =    parNN(p*d + d +(1 : d));
           xbOut =    parNN(p*d + d + d + 1);
           xbOut = repmat(xbOut,n,1);
           xhVal1 = obj.actFun(obj.X1 * xWHidden + xbHidden);
           vx = xhVal1 * xWOut + xbOut+obj.basisDVx * parHM;%+ xbOut
           
           % Calculate vy
           yWHidden = parNN(parNum + (1 : p * d));
           yWHidden = reshape(yWHidden, p, d);
           ybHidden = parNN(parNum + p*d + (1 : d));
           ybHidden = repmat(ybHidden',n,1);
           yWOut =    parNN(parNum + p*d + d +(1 : d));
           ybOut =    parNN(parNum + p*d + d + d + 1);
           ybOut = repmat(ybOut,n,1);
           yhVal1 = obj.actFun(obj.X2 * yWHidden + ybHidden);
           vy = yhVal1 * yWOut + ybOut + obj.basisDVy * parHM; %+ ybOut
           
           % Calculate vz
           zWHidden = parNN(parNum*2 + (1 : p * d));
           zWHidden = reshape(zWHidden, p, d);
           zbHidden = parNN(parNum*2 + p*d + (1 : d));
           zbHidden = repmat(zbHidden',n,1);
           zWOut =    parNN(parNum*2 + p*d + d +(1 : d));
           zbOut =    parNN(parNum*2 + p*d + d + d + 1);
           zbOut = repmat(zbOut,n,1);
           zhVal1 = obj.actFun(obj.X3 * zWHidden + zbHidden);
           vz = zhVal1 * zWOut + zbOut + obj.basisDVz * parHM; %+ zbOut
        
        end
        
		function [Q] = makeQ(obj, par)
			% MAKEQ
			%    Calculate precision matrix according to parameters
            % Get coefficients
            [vx, vy, vz] = obj.calcCoeff(par(obj.idxDVx),par(obj.idxNN));
            
			% Find DV
			DV = obj.makeDV(0);
			DVInv = obj.makeDV(1);

            % Find DRho2
			DRho2 = obj.makeDRho2(par(obj.idxRho), 0);
			DRho2Inv = obj.makeDRho2(par(obj.idxRho), 1);

			% Find DSigma
			DSigmaInv = obj.makeDSigma(par(obj.idxSigma), 1);

			% Find AH
			AH = obj.makeAH(vx, vy, vz);

			% Construct Q
			Aop = (DVInv*DRho2).^(1/2)*(DV*DRho2Inv-AH)*DSigmaInv/sqrt(4*pi);
			Q = Aop'*Aop;
            
            % Store stuff for later use in derivatives
			obj.DV = DV;
			obj.DVInv = DVInv;
			obj.DRho2 = DRho2;
			obj.DRho2Inv = DRho2Inv;
			obj.AH = AH;
			obj.DSigmaInv = DSigmaInv;
			obj.Aop = Aop;
		end

		function [AH] = makeAH(obj, vx, vy, vz)
			% MAKEAH
			%    Calculate matrix according to diffusion

			% Use mex function
			nT = size(obj.meshTV,1);
			nV = size(obj.meshVLoc,1);
			% [iVec, jVec, vVec] = SPDE.mexAH_new(nT, nV, obj.meshVLoc', obj.meshTV', obj.meshTT', vx, vy, vz);
            [iVec, jVec, vVec] = SPDE.mexAH_finalDer(nT, nV, obj.meshVLoc', obj.meshTV', obj.meshTT', vx, vy, vz, 0, vx, vy, vz, obj.startIdx, obj.rIdx, obj.rVal, length(obj.rIdx), obj.tensorPar);
			AH = sparse(iVec, jVec, vVec, nT, nT);
		end

		function [DSigma] = makeDSigma(obj, par, invert)
			% MAKEDSIGMA
			%    Calculate matrix of sigmas

			% Find current sigma function
			sigma = obj.basisSigma*par;

			% Create matrix
			dim = length(sigma);
			if(invert == 0)
				DSigma = sparse(1:dim, 1:dim, exp(sigma), dim, dim);
			else
				DSigma = sparse(1:dim, 1:dim, exp(-sigma), dim, dim);
			end
		end

		function [DRho2] = makeDRho2(obj, par, invert)
			% MAKEDRHO2
			%    Calculate matrix of rho squareds

			% Find current rho function
			rho = obj.basisRho*par;

			% Create matrix
			dim = length(rho);
			if(invert == 0)
				DRho2 = sparse(1:dim, 1:dim, exp(rho).^2, dim, dim);
			else
				DRho2 = sparse(1:dim, 1:dim, exp(-rho).^2, dim, dim);
			end
		end

		function [DV] = makeDV(obj, invert)
			% MAKEDV
			%    Make matrix with areas on diagonal
			dim = length(obj.meshArea);
			if(invert == 0)
				DV = sparse(1:dim, 1:dim, obj.meshArea, dim, dim);
			else
				DV = sparse(1:dim, 1:dim, 1./obj.meshArea, dim, dim);
			end
		end

		function [A] = triangleArea(obj, p1, p2, p3)
			% TRIANGLEAREA
			%    Calculate area of triangle defined by the three points

			% Calculate side lengths
			a = norm(p1-p2);
			b = norm(p1-p3);
			c = norm(p2-p3);

			% Calculate semiperimeter
			s = (a+b+c)/2;

			% Calculate area (overly careful about underflow)
			A = sqrt(abs(s*(s-a)*(s-b)*(s-c)));
        end
	end
end