classdef INLA
    %INLA
    % Call INLA to do different things
    
    properties
    end
    
    methods(Static)
        function [Qinv] = sparseQinv(Q)
        %sparseQinv
        % Calculate sparse inverse subset of symmetric, positive definite
        % matrix Q.
        %
        % INPUT:
        %    Q: Symmetric, positive definite matrix.
        %
        % OUTPUT:
        %    Qinv: Resulting partial inverse

        %% Generate variables that will be written to binary file
            % Get row, columns and values
            [r, c, v] = find(Q);

            % Header
            header = zeros(8,1);
            header(1) = 0;
            header(2) = size(r,1);
            header(3) = size(Q,1);
            header(4) = size(Q,2);
            header(5) = 1;
            header(6) = 1;
            header(7) = 0;
            header(8) = 1;
            h_length = 8;

        %% Write to file
            % Find a (possibly) free filename
            tmp_toInla = tempname('./');
            tmp_fromInla = strcat([tmp_toInla 'fromInla']);
            fToId = fopen(tmp_toInla, 'w');

            if(fToId == -1)
                display(sprintf('Could not open file %s',tmp_toInla));
                return
            end

            % Write length of header
            fwrite(fToId, h_length, 'int');

            % Write header
            fwrite(fToId, header, 'int');

            % Write row, col and values
            r = r-1;
            fwrite(fToId, r, 'int');
            c = c-1;
            fwrite(fToId, c, 'int');
            fwrite(fToId, v, 'double');

            % Close file to avoid problems when inla opens it
            fclose(fToId);
            
        %% Will not use costraints, but INLA will complain without it
        tmp_const = tempname('./');

        %% Run INLA to get result
            % Command
            exeName = '';
            if(ispc == 1)
                exeName = '../Shared/inla64.exe';
            end
            if(isunix == 1)
                exeName = '../Shared/inla64';
            end
            if(ismac == 1)
                exeName = '../Shared/inla.run';
            end
            com = sprintf('%s -s -m qinv %s %s %s', exeName, tmp_toInla, tmp_const, tmp_fromInla);

            % RUN
            [status, result] = system(com);

            % Check status
            if(status ~= 0)
                result
                display(sprintf('Could not run command: "%s"',com));
            end

        %% Read binary file returned
            % Open file
            fFromId = fopen(tmp_fromInla, 'r');
            if(fFromId == -1)
                delete(tmp_toInla);
                error('MyApp:INLAFailed','Could not open file %s\n',tmp_fromInla);
            end

            % Read header length
            h_length = fread(fFromId, 1, 'int');

            % Read header
            header = fread(fFromId, 8, 'int');

            % Read row, column and values
            rOut = fread(fFromId, size(r,1), 'int');
            rOut = rOut+1;
            cOut = fread(fFromId, size(c,1), 'int');
            cOut = cOut+1;
            vOut = fread(fFromId, size(v,1), 'double');

            % Close file
            fclose(fFromId);

        %% Delete both files
            % DELETE
            delete(tmp_toInla, tmp_fromInla);

        %% Assemble matrix
            Qinv = sparse(rOut, cOut, vOut);
            
        end
        
	function [Qinv] = sparseQinvStrange(Q, eZ)
        %sparseQinv
        % Calculate sparse inverse subset of symmetric, positive definite
        % matrix Q.
        %
        % INPUT:
        %    Q:  Symmetric, positive definite matrix.
	%    eZ: Add eZ rows and columns of explicit zeroes
        %
        % OUTPUT:
        %    Qinv: Resulting partial inverse

        %% Generate variables that will be written to binary file
	    % add extra zeros
	    eZ = eZ-1;
	    Z = (Q~=0);
	    Z(:, end-eZ:end) = 1;
	    Z(end-eZ:end, :) = 1;
	    Z(Q~=0) = 0;
	    Q = Q + 1e300*Z;
	    [r, c, v] = find(Q);
	    v(v==1e300) = 0;

            % Header
            header = zeros(8,1);
            header(1) = 0;
            header(2) = size(r,1);
            header(3) = size(Q,1);
            header(4) = size(Q,2);
            header(5) = 1;
            header(6) = 1;
            header(7) = 0;
            header(8) = 1;
            h_length = 8;

        %% Write to file
            % Find a (possibly) free filename
            tmp_toInla = tempname('./');
            tmp_fromInla = strcat([tmp_toInla 'fromInla']);
            fToId = fopen(tmp_toInla, 'w');

            if(fToId == -1)
                display(sprintf('Could not open file %s',tmp_toInla));
                return
            end

            % Write length of header
            fwrite(fToId, h_length, 'int');

            % Write header
            fwrite(fToId, header, 'int');

            % Write row, col and values
            r = r-1;
            fwrite(fToId, r, 'int');
            c = c-1;
            fwrite(fToId, c, 'int');
            fwrite(fToId, v, 'double');

            % Close file to avoid problems when inla opens it
            fclose(fToId);

        %% Run INLA to get result
            % Command
            exeName = '';
            if(ispc == 1)
                exeName = 'inla64.exe';
            end
            if(isunix == 1)
                exeName = 'inla64';
            end
            if(ismac == 1)
                error('Mac doesnt work');
            end
            com = sprintf('./%s -s -m qinv %s %s', exeName, tmp_toInla, tmp_fromInla);

            % RUN
            [status, result] = system(com);

            % Check status
            if(status ~= 0)
                display(sprintf('Could not run command: "%s"',com));
            end

        %% Read binary file returned
            % Open file
            fFromId = fopen(tmp_fromInla, 'r');
            if(fFromId == -1)
                delete(tmp_toInla);
                error('MyApp:INLAFailed','Could not open file %s\n',tmp_fromInla);
            end

            % Read header length
            h_length = fread(fFromId, 1, 'int');

            % Read header
            header = fread(fFromId, 8, 'int');

            % Read row, column and values
            rOut = fread(fFromId, size(r,1), 'int');
            rOut = rOut+1;
            cOut = fread(fFromId, size(c,1), 'int');
            cOut = cOut+1;
            vOut = fread(fFromId, size(v,1), 'double');

            % Close file
            fclose(fFromId);

        %% Delete both files
            % DELETE
            delete(tmp_toInla, tmp_fromInla);

        %% Assemble matrix
            Qinv = sparse(rOut, cOut, vOut);
            
        end
    end
    
end
