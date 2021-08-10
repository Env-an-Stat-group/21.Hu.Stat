classdef SphericalHarmonics < handle
	% SPHERICALHARMONICS
	%    Evaluate spherical harmonics and vector spherical harmonics
	%    at different locations
	
	methods(Static)
		function [val] = calculateHarmonicsAzEl(maxL, az, el)
			% CALCULATEHARMONICSAZEL
			%    Evaluate spherical harmonics of all orders up to l = maxL
			%    where spherical coordinates are according to Matlab's
			%    cart2sph convention
			%
			%    INPUT:
			%       maxL: {0,1,...}
			%       az: [-pi pi]
			%		el: [-pi/2 pi/2]

			% Call other function with things fixed...
			val = SPDE.SphericalHarmonics.calculateHarmonics(maxL, az, pi/2-el);
		end
		function [val, valDerPhi, valDerTheta] = calculateHarmonics(maxL, phi, theta)
			% CALCULATEHARMONICS
			%    Evaluate spherical harmonics at of all order up to l = maxL
			%    at specified azimuth and elevation coordinates
			%
			%    INPUT:
			%       maxL: {0,1,...}
			%       phi: [0 2*pi]
			%		theta: [0 pi]

			% Do some sanitation
			if(size(phi,2) > size(phi,1))
				phi = phi';
			end
			if(size(theta,2) > size(theta,2))
				theta = theta';
			end

			% Initialize storage
			val = zeros(size(phi,1), (maxL+1)^2);
			if(nargin > 1)
				valDerTheta = val;
				valDerPhi = val;
			end

			% Iterate through each order of l
			for l = 0:maxL
				% Evaluate all associated Legendre polynomials needed
				P = legendre(l, cos(theta))';

				% Find normalization factors
				N = zeros(l+1,1);
				for mI = 1:length(N)
					m = mI-1;
					N(mI) = sqrt((2*l+1)/(4*pi)*factorial(l-m)/factorial(l+m));
				end
				
				% Calculate m = 0
				val(:,(l^2+1)) = N(1)*P(:,1);

				% Calculate m ~= 0
				if(l > 0)
					% Calculate m > 0
					val(:, (l^2+2):(l^2+l+1)) = sqrt(2)*kron(N(2:end)', ones(size(phi,1),1)).*P(:,2:end).*cos(kron(1:l, phi));

					% Calculate m < 0
					val(:, (l^2+l+2):(l^2+2*l+1)) = sqrt(2)*kron(N(2:end)', ones(size(phi,1),1)).*P(:,2:end).*sin(kron(1:l, phi));
				end

				% Do derivative
				if(nargin > 1)
					if(l == 0)
						valDerTheta(:,1) = 0;
						valDerPhi(:,1) = 0;
						continue;
					end

					% Find derivative with respect to phi
					% Calculate = 0


					% Calculate m > 0
					valDerPhi(:, (l^2+2):(l^2+l+1)) = sqrt(2)*kron(N(2:end)', ones(size(phi,1),1)).*P(:,2:end).*kron(1:l, ones(size(phi,1),1)).*(-sin(kron(1:l, phi)));

					% Calculate m < 0
					valDerPhi(:, (l^2+l+2):(l^2+2*l+1)) = sqrt(2)*kron(N(2:end)', ones(size(phi,1),1)).*P(:,2:end).*sin(kron(1:l, phi));

					% Do derivative of P
					P2 = [P(:,2:end), zeros(size(P,1),1)];
					P0 = [-1/(l*(l+1))*P(:,2), P(:,1:end-1)];

					% Derivative of associated legendre part
					Pder = 0.5*P2-0.5*P0*diag((l+(0:l)).*(l-(0:l)+1));

					% Find derivative with respect to theta
					% Calculate m > 0
					valDerTheta(:, (l^2+2):(l^2+l+1)) = sqrt(2)*kron(N(2:end)', ones(size(phi,1),1)).*Pder(:,2:end).*cos(kron(1:l, phi));

					% Calculate m < 0
					valDerTheta(:, (l^2+l+2):(l^2+2*l+1)) = sqrt(2)*kron(N(2:end)', ones(size(phi,1),1)).*Pder(:,2:end).*sin(kron(1:l, phi));
				end
			end
		end

		function [Yhar, vx, vy, vz] = sphericalR3(l, phi, theta)
			% Get stuff in \phi, \theta coordinates
			[Yhar, YderPhi, YderTheta] = SPDE.SphericalHarmonics.sphericalReal(l, phi, theta);

			% Initialize storage
			vx = zeros(size(YderPhi, 1), 2*size(YderPhi, 2));
			vy = zeros(size(YderPhi, 1), 2*size(YderPhi, 2));
			vz = zeros(size(YderPhi, 1), 2*size(YderPhi, 2));

			% Vector field 1
			vecTheta = YderTheta;
			vecPhi = YderPhi;

			% Convert to x, y, z
			len = size(YderPhi,2);
			vx(:, 1:len) = kron(ones(1, len), cos(theta).*cos(phi)).*vecTheta + kron(ones(1, len), -sin(phi)).*vecPhi;
			vy(:, 1:len) = kron(ones(1, len), cos(theta).*sin(phi)).*vecTheta + kron(ones(1, len), cos(phi)).*vecPhi;
			vz(:, 1:len) = kron(ones(1, len), -sin(theta)).*vecTheta;

			% Vector field 2
			vecTheta = -YderPhi;
			vecPhi = YderTheta;

			% Convert to x, y, z
			vx(:, (len+1):end) = kron(ones(1, len), cos(theta).*cos(phi)).*vecTheta + kron(ones(1, len), -sin(phi)).*vecPhi;
			vy(:, (len+1):end) = kron(ones(1, len), cos(theta).*sin(phi)).*vecTheta + kron(ones(1, len), cos(phi)).*vecPhi;
			vz(:, (len+1):end) = kron(ones(1, len), -sin(theta)).*vecTheta;
		end

		function [Yhar, YderPhi, YderTheta] = sphericalReal(l, phi, theta)
			% Create index vector
			m = 1:l;

			% Get complex spherical harmonics
			YharTmp = SPDE.SphericalHarmonics.sphericalComplex(l, phi, theta);

			% Make real
			Yhar = YharTmp;
			negM = fliplr(m);
			posM = l+1+m;
			Yhar(:,negM) = 1/(1.0i*sqrt(2))*(YharTmp(:, negM)-YharTmp(:, posM)*diag((-1).^m));			
			Yhar(:,posM) = 1/sqrt(2)*(YharTmp(:,posM)+YharTmp(:,negM)*diag((-1).^m));

			% Make certain it is real
			Yhar = real(Yhar);

			% Maybe do derivatives
			if(nargout > 1)
				% Evaluate again...
				[~, YderPhiTmp, YderThetaTmp] = SPDE.SphericalHarmonics.sphericalComplex(l, phi, theta);

				% Fix derivative w.r.t. phi
				YderPhi = YderPhiTmp;
				YderPhi(:,negM) = 1/(1.0i*sqrt(2))*(YderPhiTmp(:, negM)-YderPhiTmp(:, posM)*diag((-1).^m));			
				YderPhi(:,posM) = 1/sqrt(2)*(YderPhiTmp(:,posM)+YderPhiTmp(:,negM)*diag((-1).^m));

				% Remove complex part (should be 0)
				YderPhi = real(YderPhi);

				% Fix derivative w.r.t. theta
				YderTheta = YderThetaTmp;
				YderTheta(:,negM) = 1/(1.0i*sqrt(2))*(YderThetaTmp(:, negM)-YderThetaTmp(:, posM)*diag((-1).^m));			
				YderTheta(:,posM) = 1/sqrt(2)*(YderThetaTmp(:,posM)+YderThetaTmp(:,negM)*diag((-1).^m));

				% Remove complex part(should be 0)
				YderTheta = real(YderTheta);
			end
		end

		function [Yhar, YderPhi, YderTheta] = sphericalComplex(l, phi, theta)
			% m-values
			m = (-l):l;

			% Normalizing constant
			C = sqrt((2*l+1)/(4*pi)*factorial(l-m)./factorial(l+m));

			% Calculate associated legendre
			P = legendre(l, cos(theta))';

			% Find negative ones
			mpos = 1:l;
			C2 = ((-1).^mpos).*factorial(l-mpos)./factorial(l+mpos);
			P = [fliplr(P(:,2:end)*diag(C2)), P];

			% Calculate
			Yhar = P.*exp(1.0i*kron(m,phi))*diag(C);

			if(nargout > 1)
				% Calculate derivative w.r.t. phi
				lNext = l+1;
				mNext = (-lNext):lNext;
				CNext = sqrt((2*lNext+1)/(4*pi)*factorial(lNext-mNext)./factorial(lNext+mNext));
				PNext = legendre(lNext, cos(theta))';
				mposNext = 1:lNext;
				C2Next = ((-1).^mposNext).*factorial(lNext-mposNext)./factorial(lNext+mposNext);
				PNext = [fliplr(PNext(:,2:end)*diag(C2Next)), PNext];
				
				% Do recursion to get rid of 1/sin(theta)
				YderPhi = -0.5*(PNext(:,3:end)+PNext(:,1:end-2)*diag((l-m+1).*(l-m+2)))*diag(1./m);
				YderPhi(:,l+1) = 0;
				YderPhi = YderPhi.*exp(1.0i*kron(m,phi))*diag(C)*diag(1.0i*m);

				% Pre-do som stuff
				P2 = [P(:,2:end), zeros(size(P,1),1)];
				P0 = [zeros(size(P,1),1), P(:,1:end-1)];

				% Find derivative of P
				Pder = 0.5*P2-0.5*P0*diag((l+m).*(l-m+1));

				% Calc derivative w.r.t. theta
				YderTheta = Pder.*exp(1.0i*kron(m,phi))*diag(C);
			end
		end
	end
end
