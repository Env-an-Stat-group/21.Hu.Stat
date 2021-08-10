classdef MeshGenerator < handle
	% MESHGENERATOR
	%    Create triangulation of a 2D area or a sphere

	properties
	end

	methods(Static)
		function [vLoc, tt, tv] = generateMesh(loc, varargin)
			% GENERATEMESH
			%    Generate the triangular mesh

			% Write FMesher file
			folder = 'DontPutFilesInThisFolder';
			prefix = strcat([folder '/mesh.']);
			fname = 'input.s';
			totfname = strcat([prefix fname]);
			mkdir(folder);
			SPDE.MeshGenerator.writeFMesherFile(loc, totfname);

			% Set some default values that doesn't give insane grid sizes
			Ar = prod(range(loc));
			e1 = 16;
			e2 = -0.1;
			angl = 21;
			desNumTriangle = 10000;
			E1 = sqrt(Ar/(sqrt(3)/4*desNumTriangle))*sqrt(2.5);
			E2 = E1;
			CO = E1;
			globe = 0;

			% Parse optional input
			p = inputParser;
			addOptional(p, 'extend1', e1);
			addOptional(p, 'extend2', e2);
			addOptional(p, 'minAngle', angl);
			addOptional(p, 'maxEdge1', E1);
			addOptional(p, 'maxEdge2', E2);
			addOptional(p, 'cutOff', CO);
			addOptional(p, 'numTriangles', desNumTriangle);
			addOptional(p, 'globe', globe);
			parse(p, varargin{:});
			e1 = p.Results.extend1;
			e2 = p.Results.extend2;
			angl = p.Results.minAngle;
			E1 = p.Results.maxEdge1;
			E2 = p.Results.maxEdge2;
			CO = p.Results.cutOff;
			desNumTriangle = p.Results.numTriangles;
			globe = p.Results.globe;
			if(desNumTriangle ~= 10000)
				E1 = sqrt(Ar/(sqrt(3)/4*desNumTriangle))*sqrt(2.5);
				E2 = E1;
				CO = E1;
			end

			% Set parameter values
			extend1 = sprintf('%f', e1);
			extend2 = sprintf('%f', e2);
			cutOff = sprintf('%f', CO);
			ang = sprintf('%f', angl);
			maxE1 = sprintf('%f', E1);
			maxE2 = sprintf('%f', E2);

			% Check operating system
			exeName = '';
			if(isunix == 1)
				exeName = '../Shared/fmesher64';
			end
			if(ispc == 1)
				exeName = '../Shared/fmesher64.exe';
			end
			if(ismac == 1)
				exeName = '../Shared/fmesher.run';
			end

			% Check if triangulation of 2-sphere should be generated
			if(globe > 0)
				globe = sprintf('%d', globe);
				cmd = strcat([exeName ' --globe=' globe ' --cutoff=' cutOff ' --rcdt=' ang ',' maxE1 ',' maxE2 ' ' prefix]);
			else
				cmd = strcat([exeName ' --input=' fname ' --cutoff=' cutOff ' --cet=' extend1 ',' extend2 ' --rcdt=' ang ',' maxE1 ',' maxE2 ' ' prefix]);
			end	

			% Run command
		    display(cmd)
			system(cmd);

			% Get vertex locations
			fname = 's';
			totfname = strcat([prefix fname]);
			vLoc = SPDE.MeshGenerator.readFMesherFile(totfname);

			% Get triangle-triangle map
			fname = 'tt';
			totfname = strcat([prefix fname]);
			tt = SPDE.MeshGenerator.readFMesherFile(totfname);
			tt(tt>=0) = tt(tt>=0)+1;

			% Get triangle-vertex map
			fname = 'tv';
			totfname = strcat([prefix fname]);
			tv = SPDE.MeshGenerator.readFMesherFile(totfname);
			tv = tv+1;

			% Remove temporary directory
			%rmdir(folder, 's');
		end


		function [] = writeFMesherFile(loc, fname)
			% WRITEFMESHERFILE
			%    Write a file to be used by fmesher binary

			% Extract required information
			nrow = size(loc, 1);
			ncol = size(loc, 2);
			elems = nrow*ncol;
			datatype = 0;
			valuetp = 1;
			matrixtype = 0;
			storagetype=1;
			vers = 0;

			% Collect in header vector
			h = [vers; elems; nrow; ncol; datatype; valuetp; matrixtype; storagetype];

			% Open file
			fp = fopen(fname, 'w');

			% Write header size
			fwrite(fp, length(h), 'int');

			% Write header
			fwrite(fp, h, 'int');

			% Write locations
			fwrite(fp, loc(:), 'double');

			% Close file
			fclose(fp);
		end


		function [val] = readFMesherFile(fname)
			% READFMESHERFILE
			%    Read a file returned by fmesher binary

			% Open file
			fp = fopen(fname);

			% read header length (and ignore...)
			fread(fp, 1, 'int');

			% Read header
			h = fread(fp, 8, 'int');

			% extract interesting information
			nrow = h(3);
			ncol = h(4);
			valuetp = h(6);
			storagetype = h(8);

			% Read int/double(s)
			if(valuetp == 0)
				val = fread(fp, nrow*ncol, 'int');
			else
				val = fread(fp, nrow*ncol, 'double');
			end
			
			% Adjust according to storage type
			if(storagetype == 0)
				val = reshape(val, [ncol nrow])';
			else
				val = reshape(val, [nrow ncol]);
			end
		end
	end
end	
