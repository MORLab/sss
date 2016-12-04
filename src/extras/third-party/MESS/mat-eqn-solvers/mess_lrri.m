function [ Z, out ] = mess_lrri( eqn, opts, oper )
% MESS_LRRI
% Solve the stable generalized algebraic H-infinity Riccati equation
%
%   (N) BB' + AXE' + EXA' + EX(C1'C1 - C2'C2)XE' = 0
%   (T) C'C + A'XE + E'XA + E'X(B1B1' - B2B2')XE = 0
%
% for a full-rank factor of X via Riccati iteration.
% Using directly and indirectly Matlab M.E.S.S. structures and functions
% for solving definite Riccati equations.
%
% INPUT:
%   eqn  - a data structure containing the data of the H-infinity Riccati 
%          equation to solve, with following possible entries:
%          (for example eqn.A_)
%   |---------|-----------------------------------------------------------|
%   |  ENTRY  |                         MEANING                           |
%   |---------|-----------------------------------------------------------|
%   | A_      | Matrix or Operator from (N), (T) corresponding to the     |
%   |         | user defined operations like mul_A() loaded in the oper   |
%   |         | structure                                                 |
%   |---------|-----------------------------------------------------------|
%   | B       | a (n x m)-matrix from (N), only used if eqn.type == 'N'   |
%   |---------|-----------------------------------------------------------|
%   | B1      | a (n x m1)-matrix from (T), containing the disturbances,  |
%   |         | only used if eqn.type == 'T'                              |
%   |---------|-----------------------------------------------------------|
%   | B2      | a (n x m2)-matrix from (T), containing the controls, only |
%   |         | used if eqn.type == 'T'                                   |
%   |---------|-----------------------------------------------------------|
%   | C       | a (p x n)-matrix from (T), only used if eqn.type == 'T'   |
%   |---------|-----------------------------------------------------------|
%   | C1      | a (p1 x n)-matrix from (N), containing the errors, only   |
%   |         | used if eqn.type == 'N'                                   |
%   |---------|-----------------------------------------------------------|
%   | C2      | a (p2 x n)-matrix from (N), containing the outputs, only  |
%   |         | used if eqn.type == 'N'                                   |
%   |---------|-----------------------------------------------------------|
%   | E_      | Matrix or Operator from (N), (T) corresponding to the     |
%   |         | user defined operations like mul_E() loaded in the oper   |
%   |         | structure                                                 |
%   |---------|-----------------------------------------------------------|
%   | type    | containg either 'N' to solve equation (N) or otherwise    |
%   |         | 'T' to solve equation (T)                                 |
%   |---------|-----------------------------------------------------------|
%
%   opts - an option structure containing all desired options for each used
%          method, i.e., opts.adi contains options for mess_lradi, opts.nm
%          for mess_lrnm and opts.ri for mess_lrri. Here only the structure
%          opts.ri is considered with following entries:
%          (for example opts.ri.maxiter)
%   |---------|-----------------------------------------------------------|
%   |  ENTRY  |                         MEANING                           |
%   |---------|-----------------------------------------------------------|
%   | maxiter | nonnegative scalar, upper bound for iteration steps of    |
%   |         | Riccati iteration                                         |
%   |---- ----|-----------------------------------------------------------|
%   | colctol | nonnegative scalar, relative accuracy for column          |
%   |         | compression with singular value decomposition             |
%   |---------|-----------------------------------------------------------|
%   | info    | { 0, 1 } for activating display mode                      |
%   |---------|-----------------------------------------------------------|
%   | restol  | nonnegative scalar, upper bound for normalized residuum,  |
%   |         | used as convergence criterion                             |
%   |---------|-----------------------------------------------------------|
%
%   oper - a function structure containing all necessary operations with
%          eqn.A_ and eqn.E_, result of the function operatormanager()
%          WARNING: tested only with 'default' usfs
%
% OUTPUT:
%   Z   - numerical full-rank factor of X from (N) or (T), such that X=Z*Z'
%   out - an information structure containing informations about all used
%         methods. The last use of mess_lradi is saved in out.adi, the last
%         use of mess_lrnm is saved in out.nm and all informations about
%         mess_lrri are saved in out.ri. Here only the structure out.ri is 
%         considered with following entries:
%         ( for example out.ri.niter)
%   |---------|-----------------------------------------------------------|
%   |  ENTRY  |                         MEANING                           |
%   |---------|-----------------------------------------------------------|
%   | niter   | number of performed Riccati iteration steps               |
%   |---------|-----------------------------------------------------------|
%   | res     | vector, containing the computed normalized rediduals for  |
%   |         | every performed iteration step                            |
%   |---------|-----------------------------------------------------------|
%
%
% NOTE: The operator functions in usfs, result of the function 
%       operatormanager(), are tested only for the 'default' case. Usage of
%       others at own risk.
%
%
% Steffen Werner, 2015-10-12.

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
%

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check operatorfunctions in oper                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (~isfield(oper, 'name')) || (~strcmp(oper.name, 'default'))
        warning( 'MESS:untested', ['The operator functions in usfs, ',...
            'result of the function operatormanager(), are tested only for' ,...
            ' the ''default'' case. Usage of others at own risk.'] );
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for ADI Control structure in options                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Single fields are checked below or inside mess_lradi.
    if (~isfield( opts, 'adi' )) || (~isstruct( opts.adi ))
        error( 'MESS:control_data', 'ADI control structure opts.ADI missing.' );
    end
    
    if ~isfield( opts.adi, 'computeZ' )
        opts.adi.computeZ = 1; 
    end
    
    if ~isfield( opts.adi, 'accumulateK' )
        opts.adi.accumulateK = 0; 
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for shift parameter structure                                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (~isfield( opts.adi, 'shifts' )) || (~isstruct(opts.adi.shifts))
        error( 'MESS:control_data', 'shift parameter control structure missing.' );
    else
        if (~isfield( opts.adi.shifts, 'l0' )) || (~isfield( opts.adi.shifts, 'kp' )) || (~isfield( opts.adi.shifts, 'km' ))
            error( 'MESS:shifts', 'Incomplete input parameters for shift computation.' );
        end  
    end
    
    if ~isfield( opts.adi.shifts, 'period' )
        opts.adi.shifts.period = 1;
    end
    
    if ~isfield( opts.adi.shifts, 'method' )
        opts.adi.shifts.method = 'heur';
    end
    
    if strcmp( opts.adi.shifts.method, 'wachspress' )
        if ~isfield( opts.adi.shifts, 'wachspress' )
            opts.adi.shifts.wachspress = 'T';
        end
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for Newton control structure in options                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (~isfield( opts, 'nm' )) || (~isstruct( opts.nm ))
        error( 'MESS:control_data', 'Newton control structure opts.nm missing.' );
    else
        
        if (~isfield( opts.nm, 'maxiter' )) || (~isnumeric( opts.nm.maxiter ))
            warning( 'MESS:maxiter', 'Missing or Corrupted maxiter field. Switching to default.' );
            opts.nm.maxiter = 20;
        end

        if (~isfield( opts.nm, 'rctol' )) || (~isnumeric( opts.nm.rctol ))
            warning( 'MESS:rctol', 'Missing or Corrupted rctol field. Switching to default.' );
            opts.nm.rctol = 0;
        end

        if (~isfield( opts.nm, 'restol' )) || (~isnumeric( opts.nm.restol ))
            warning( 'MESS:restol', 'Missing or Corrupted restol field. Switching to default.' );
            opts.nm.restol = 0;
        end

        if (~isfield( opts.nm, 'info' )) || (~isnumeric( opts.nm.info ))
            warning( 'MESS:info', 'Missing or Corrupted info field. Switching to default.' );
            opts.nm.info = 0;
        end
      
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for Riccati control structure in options                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (~isfield( opts, 'ri' )) || (~isstruct( opts.ri ))
        error( 'MESS:control_data', 'Riccati control structure opts.ri missing.' );
    else

        if (~isfield( opts.ri, 'maxiter' )) || (~isnumeric( opts.ri.maxiter ))
            warning( 'MESS:maxiter', 'Missing or Corrupted maxiter field. Switching to default.' );
            opts.ri.maxiter = 10;
        end
        
        if (~isfield( opts.ri, 'restol' )) || (~isnumeric( opts.ri.restol ))
            warning( 'MESS:restol', 'Missing or Corrupted restol field. Switching to default.' );
            opts.ri.restol = 0;
        end
        
        if (~isfield( opts.ri, 'colctol' )) || (~isnumeric( opts.ri.colctol ))
            warning( 'MESS:colctol', 'Missing or Corrupted colc field. Switching to default.' );
            opts.ri.colctol = 0;
        end
        
        if (~isfield( opts.ri, 'info' )) || (~isnumeric( opts.ri.info ))
            warning( 'MESS:info', 'Missing or Corrupted info field. Switching to default.' );
            opts.ri.info = 0;
        end
        
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check system data                                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield( eqn, 'type' )
        eqn.type = 'N';
        warning( 'MESS:equation_type', 'Unable to determine type of equation. Falling back to type ''N''.' );
    elseif (eqn.type ~= 'N') && (eqn.type ~= 'T')
        error( 'MESS:equation_type', 'Equation type must be either ''T'' or ''N''.' );
    end
    
    if ~isfield( eqn, 'haveE' )
        eqn.haveE = 0; 
    end
    
    [ eqn, erg ] = oper.init( eqn, opts, 'A', 'E' );
    
    if ~erg
        error( 'MESS:control_data', 'system data is not completely defined or corrupted' );
    end
    
    % make sure the corresponding matrices from quadratic term are well
    % defined and the first right hand side is dense so that the resulting 
    % factor is densly stored.
    if eqn.type == 'T'
        if (~isfield( eqn, 'B1' )) || (~isnumeric( eqn.B1 ))
            error( 'MESS:control_data', 'eqn.B1 is not defined or corrupted' );
        end
        
        if (~isfield( eqn, 'B2' )) || (~isnumeric( eqn.B2 ))
            error( 'MESS:control_data', 'eqn.B2 is not defined or corrupted' );
        end
        
        if (~isfield( eqn, 'C' )) || (~isnumeric( eqn.C ))
            error( 'MESS:control_data', 'eqn.C is not defined or corrupted' );
        end
        
        if issparse( eqn.B1 )
            eqn.B1 = full( eqn.B1 ); 
        end
        
        if issparse( eqn.B2 )
            eqn.B2 = full( eqn.B2 ); 
        end
        
        if issparse( eqn.C )
            eqn.C = full( eqn.C ); 
        end
    else
        if (~isfield( eqn, 'C1' )) || (~isnumeric( eqn.C1 ))
            error( 'MESS:control_data', 'eqn.C1 is not defined or corrupted' );
        end
        
        if (~isfield( eqn, 'C2' )) || (~isnumeric( eqn.C2 ))
            error( 'MESS:control_data', 'eqn.C2 is not defined or corrupted' );
        end
        
        if (~isfield( eqn, 'B' )) || (~isnumeric( eqn.B ))
            error( 'MESS:control_data', 'eqn.B is not defined or corrupted' );
        end
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                     %
    % All checks done. Here comes the real work!                          %
    %                                                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize data                                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.ri.restol
        res = zeros( opts.ri.maxiter, 1 ); 
    else
        res = [];
    end

    if eqn.type == 'T'
        eqn.B = eqn.B2;
    else
        eqn.C = eqn.C2;
    end
    
    if opts.ri.restol
        if eqn.type == 'T'
            res0  = norm( eqn.C * eqn.C', 2 );
        else
            res0  = norm( eqn.B' * eqn.B, 2 );
        end
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize usfs                                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ eqn, opts, oper ] = oper.mul_E_pre( eqn, opts, oper );
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start iteration                                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = 1;
    Z = [];
    
    while j <= opts.ri.maxiter
        
        % Solve updated Riccati equation for full-rank factor.
        [ Y, nmout ] = mess_lrnm( eqn, opts, oper );
        
        % Remove initial solution from equation.
        if isfield( eqn, 'K0' )
            eqn = rmfield( eqn, 'K0' );
        end
        
        % Column compression via SVD.
        [ U, S, ~ ] = svd( [ Z, Y ], 'econ' );
        r           = find( diag( S ) > S(1,1) * opts.ri.colctol, 1, 'last' );
        Z           = U * S(:,1:r);
        
        % Set rank-k update parameter.
        eqn.setUV = 1;
        
        % Update the constant term and error variables.
        if eqn.type == 'T'
            if eqn.haveE
                eqn.C = (oper.mul_E( eqn, opts, 'T', Y, 'N' ) * (Y' * eqn.B1))';
            else
                eqn.C = (eqn.B1' * Y) * Y';
            end
            
            if opts.ri.restol
                res(j) = norm( eqn.C * eqn.C', 2 ) / res0;
            end
        else
            if eqn.haveE
                eqn.B = oper.mul_E( eqn, opts, 'T', Y, 'N' ) * (eqn.C1 * Y)';
            else
                eqn.B = Y * (eqn.C1 * Y)';
            end
            
            if opts.ri.restol
                res(j) = norm( eqn.B' * eqn.B, 2 ) / res0;
            end
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % print status information                                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if opts.ri.info
            if opts.nm.restol
                fprintf( 1, [ '\nRI step: %4d  normalized residual: %e\n' ...
                              '               number of NM steps: %d\n' ...
                              '               number of ADI steps: %d \n\n' ], ...
                         j, res(j), nmout.nm.niter, nmout.adi.niter );
            end
        end
        
        % Set rank-k update.
        if eqn.type == 'T'
            if eqn.haveE
                eqn.V = oper.mul_E( eqn, opts, 'T', Z, 'N' );
            else
                eqn.V = Z;
            end
            eqn.U      = eqn.B2 * (eqn.B2' * Z) - eqn.B1 * (eqn.B1' * Z);
            eqn.UVsize = size( eqn.V, 2 );
        else
            if eqn.haveE
                eqn.U = oper.mul_E( eqn, opts, 'N', Z, 'N' );
            else
                eqn.U = Z;
            end
            eqn.V      = eqn.C2' * (eqn.C2 * Z) - eqn.C1' * (eqn.C1 * Z);
            eqn.UVsize = size( eqn.U, 2 );
        end
        
        % Update iteration variable.
        j = j + 1;
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate stopping criteria                                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (opts.ri.restol) && (res(j-1) < opts.ri.restol)
            break;
        end     
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % finalize usfs                                                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ ~, opts, ~ ] = oper.mul_E_post( eqn, opts, oper );
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prepare output                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out.adi      = nmout.adi;
    out.nm       = nmout.nm;
    out.ri.niter = j - 1;
    
    if opts.ri.restol
        out.ri.res = res(1:out.ri.niter); 
    end

end
