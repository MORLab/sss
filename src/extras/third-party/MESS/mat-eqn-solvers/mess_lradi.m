function [Z,out,eqn,opts,oper]=mess_lradi(eqn,opts,oper)
%% function [Z,out,eqn,opts,oper]=mess_lradi(eqn,opts,oper)
%
% Solve continuous-time Lyapunov equations with sparse coefficients
%   eqn.type = 'N' -> A*Z*Z'*E' + E*Z*Z'*A' + B*B' = 0 (N)
%   eqn.type = 'T' -> A'*Z*Z'*E + E'*Z*Z'*A + C'*C = 0 (T)
%   
%
% Input
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation 
%                       with A and E
%
% Output
%   Z                   low rank solution factor
%
%   out                 struct containing output information
%
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation 
%                       with A and E
%
% Input fields in struct eqn:
%   eqn.A_      sparse (n x n) matrix A
%
%   eqn.E_      sparse (n x n) matrix E
%
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
%
%   eqn.G       dense (n x m1) matrix G
%               if present it is used instead of B as RHS
%               (required for LDL^T formulation otherwise optional)
%
%   eqn.S       dense (m1 x m1) matrix 
%               (required for LDL^T formulation)
%
%   eqn.U       dense (n x m3) matrix U
%               (required if eqn.V is present)
%
%   eqn.V       dense (n x m3) matrix V
%               (required if eqn.U is present)
%
%   eqn.type    possible  values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional)
%
%   eqn.haveE   possible  values: 0, 1, false, true
%               if haveE = 0: matrix E in eqn.E_ is assumed to be identity
%               (optional)
%
%   eqn.setUV   possible  values: 0, 1, false, true
%               if setUV = 1: U = [U1, U2] and V = [V1, V2]
%               if K or DeltaK are accumulated during the iteration they
%               use only U1 and V1
%               (optional)
%
% Input fields in struct opts:
%   opts.adi.maxiter            possible  values: integer > 0
%                               maximum iteration number
%                               (optional)
%
%   opts.adi.restol             possible  values: scalar >= 0
%                               stopping tolerance for the relative
%                               residual norm; if restol = 0 the relative
%                               residual norm is not evaluated
%                               (optional)
%
%   opts.adi.rctol              possible  values: scalar >= 0
%                               stopping tolerance for the relative
%                               change of the solutiuon Z; 
%                               if restol = 0 the relative
%                               change is not evaluated
%                               (optional)
%
%   opts.adi.info               possible  values: 0, 1, false, true
%                               turn on (1) or off (0) the status output in
%                               every iteration step
%                               (optional)
%
%   opts.adi.norm               possible  values: 2, 'fro'
%                               use 2-norm (2) or Frobenius norm ('fro') to
%                               compute residual and relative change norms
%                               (optional)
%
%   opts.adi.computeZ           possible  values: 0, 1, false, true
%                               turn on (1) or off (0) the computation of
%                               the solution Z; turm off if only the 
%                               feedbackmatrix K is of interest
%                               (optional)
%
%   opts.adi.accumulateK        possible  values: 0, 1, false, true
%                               accumulate the feedback matrix K during the
%                               iteration
%                               (optional)
%
%   opts.adi.accumulateDeltaK   possible  values: 0, 1, false, true
%                               accumulate the update DeltaK of the 
%                               feedback matrix K during the iteration
%                               (optional)
%
%   opts.adi.LDL_T              possible  values: 0, 1, false, true
%                               use LDL^T formulation for the RHS and
%                               solution
%                               (optional)
%
%   opts.adi.shifts.l0          possible  values: integer > 0
%                               number of shifts that should be computed
%                               2*l0 < kp + km is required
%                               (optional)
%
%   opts.adi.shifts.kp          possible  values: integer > 0
%                               number of Arnoldi steps w.r.t. A for
%                               heuristic shift computation
%                               kp < n is required
%                               (optional)
%
%   opts.adi.shifts.km          possible  values: integer > 0
%                               number of Arnoldi steps w.r.t. inv(A) for
%                               heuristic shift computation
%                               km < n is required
%                               (optional)
%
%   opts.adi.shifts.b0          (n x 1) array
%                               start vector for Arnoldi algorithm for
%                               heuristic shift computation
%                               (optional)
%
%   opts.adi.shifts.p           array with shifts
%                               complex shifts are possible
%
%   opts.adi.shifts.info        possible  values: 0, 1, false, true
%                               turn output of used shifts before the first
%                               iteration step on (1) or off (0) 
%                               (optional)
%
%   opts.adi.shifts.method      possible  values: 'heur','heuristic',
%                               'penzl','Penzl', 'wachspress','Wachspress',
%                               'projection'
%                               method for shift computation
%                               in case of 'projection' new shifts are
%                               computed during the iteration steps
%                               (optional)
%
%   opts.adi.projection.freq    possible  values: integer > 0
%                               frequency of the usage of galerkin
%                               projection acceleration
%                               (optional)
%
%   opts.adi.projection.ortho   possible  values: 0, 1, false, true
%                               implicit (0) or explicit (1) 
%                               orthogonalization of Z
%                               (optional)
%
%   opts.adi.projection.meth    possible  values: 'lyapchol', 
%                               'lyap_sgn_fac', 'lyap', 'lyapunov',
%                               'lyap2solve'
%                               projection method
%                               (optional)
%                               
%
% If optional input arguments are missing they may be set to default values
% and a 'MESS:control_data' warning is printed. to turn warnings off use
% warning('OFF', 'MESS:control_data')
%
% Matrix A can have the form A = Ã + U*V'
%     if U (eqn.U) and V (eqn.V) are provided
%     U and V are dense (n x m3) matrices and shoud have low rank m3 << n
%
% The feedback matrix K can be accumulated during the iteration:
%     eqn.type = 'N' -> K = E*Z*Z'*C'
%     eqn.type = 'T' -> K = E'*Z*Z'*B
%
% For LDL^T formulation use opts.adi.LDL_T = 1:
%     A*L*D*L'*E' + E*L*D*L'*A' + G*S*G' = 0
%     RHS has form G * S * G'
%     Solution has form L * D * L' with D Kronecker product of out.D and S
%     D is not computed explicitly
%     L is stored in Z if computed (opts.adi.computeZ)
%     G (eqn.G) and S (eqn.S) need to be given
%
% Output fields in struct out:
%   out.D               scalar factor vector for LDL^T formulation
%                       use kron(diag(out.D),eqn.S) to build D
%                       (opts.adi.LDL_T = 1)
%
%   out.res             array of relative residual norms
%
%   out.rc              array of relative change norms
%
%   out.niter           number of ADI iterations
%
%   out.res_fact        low rank residual factor W
%
%   out.Riccati_res     outer Riccati residual norm for Newton iteration
%                       (opts.nm.accumulateRes = 1)
%
%   out.linesearch      flag to trigger linesearch in Newton iteration
%                       (opts.adi.inexact ~= 0)
%
%   out.restart         flag to trigger complete restart of Newton
%                       iteration because of divergence
%
% uses oparatorfunctions size, init, init_res, init_res_pre, init_res_post,
% init_res_post, sol_ApE, mul_E, mul_E_pre, mul_E_post

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
% Copyright (C) Jens Saak, Martin Koehler and others 
%               2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
%

%% check field opts.adi
if ~isfield(opts,'adi') || ~isstruct(opts.adi)
    error('MESS:control_data',['No adi control data found in options', ...
        'structure.']);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for shifts and their properness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
illelgal_shifts=0;
if ~isfield(opts.adi,'shifts')||~isstruct(opts.adi.shifts)
    error('MESS:control_data','shift parameter control structure missing.');
end
if ~(isfield(opts.adi.shifts,'p')&&isnumeric(opts.adi.shifts.p)&&isvector(opts.adi.shifts.p))
    error('MESS:shifts', ...
        'Found empty shift vector. Please provide proper shifts.');
else
    % Check if all shifts are in the open left half plane
    if any(~(real(opts.adi.shifts.p)<0)), illelgal_shifts=1; end
    
    % Check if complex pairs of shifts are properly ordered.
    i=1;
    while i<=length(opts.adi.shifts.p)
        if ~(isreal(opts.adi.shifts.p(i)))
            if (opts.adi.shifts.p(i+1)~=conj(opts.adi.shifts.p(i))), illelgal_shifts=1;end
            i=i+1;
        end
        i=i+1;
    end
end
if illelgal_shifts
    error('MESS:shifts_improper','Improper shift vector detected!')
end
if isfield(opts.adi.shifts, 'method') && ...
        strcmp(opts.adi.shifts.method, 'projection')
    opts.adi.computeZ = 1;
    opts.adi.shifts.used_shifts = [];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check info parameter for output verbosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts.adi,'info')
    opts.adi.info=0;
else
    if ~isnumeric(opts.adi.info)&&~islogical(opts.adi.info)
        error('MESS:info','opts.adi.info parameter must be logical or numeric.')
    end
end

if ~isfield(opts.adi.shifts,'info')
    opts.adi.shifts.info=0;
else
    if ~isnumeric(opts.adi.shifts.info)&&~islogical(opts.adi.shifts.info)
        error('MESS:info','opts.adi.shifts.info parameter must be logical or numeric.')
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check stopping parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts.adi,'maxiter')||~isnumeric(opts.adi.maxiter)
    warning('MESS:control_data',...
        ['Missing or Corrupted opts.adi.maxiter field.', ...
        'Switching to default: 100']);
    opts.adi.maxiter=100;
end

if ~isfield(opts.adi,'rctol')||~isnumeric(opts.adi.rctol)
    warning('MESS:control_data',...
        ['Missing or Corrupted opts.adi.rctol field.', ...
        'Switching to default: 0']);
    opts.adi.rctol=0;
end
if opts.adi.rctol
    nrmZ = 0;
end

if ~isfield(opts.adi,'restol')||~isnumeric(opts.adi.restol)
    warning('MESS:control_data',...
        ['Missing or Corrupted opts.adi.restol field.', ...
        'Switching to default: 0']);
    opts.adi.restol=0;
end
if ~isfield(opts.adi, 'norm') || ...
        (~strcmp(opts.adi.norm, 'fro') && ...
        (~isnumeric(opts.adi.norm) || opts.adi.norm ~= 2))
    warning('MESS:control_data', ...
        ['Missing or Corrupted opts.adi.norm field.', ...
        'Switching to default: ''fro''']);
    opts.adi.norm = 'fro';
end
if ~isfield(opts.adi,'inexact'),  opts.adi.inexact=0; end
if opts.adi.inexact
    if ~opts.adi.restol
        % restol is needed
        opts.adi.restol = 1e-16;
        opts.adi.accumulateDeltaK = 1;
    end
    if ~isfield(opts.adi,'outer_tol')
        error('MESS:outer_tol',...
            'For inexact ADI opts.adi.outer_tol is needed.');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(eqn,'type')
    eqn.type='N';
    warning('MESS:control_data',['Unable to determine type of equation.'...
        'Falling back to type ''N''']);
elseif (eqn.type~='N')&&(eqn.type~='T')
    error('MESS:equation_type','Equation type must be either ''T'' or ''N''');
end

%set flag 0 if E does not exist
if ~isfield(eqn,'haveE')
    eqn.haveE=0; 
    warning('MESS:control_data', ...
        ['Missing or Corrupted eqn.haveE field.', ...
        'Switching to default: 0']);
end

[eqn, erg] = oper.init(eqn, opts, 'A','E');
if ~erg
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end

if eqn.type == 'N' && (isfield(opts.adi, 'accumulateDeltaK') && opts.adi.accumulateDeltaK)
    if ~isfield(eqn, 'B') || ~isnumeric(eqn.B)
        error('MESS:control_data', 'eqn.B is not defined or corrupted');
    end
    if ~isfield(eqn, 'C') || ~isnumeric(eqn.C)
        error('MESS:control_data', 'eqn.C is not defined or corrupted');
    end
    m = size(eqn.C, 1);
end

if eqn.type == 'T' && (isfield(opts.adi, 'accumulateDeltaK') && opts.adi.accumulateDeltaK)
    if ~isfield(eqn, 'C') || ~isnumeric(eqn.C)
        error('MESS:control_data', 'eqn.C is not defined or corrupted');
    end
    if ~isfield(eqn, 'B') || ~isnumeric(eqn.B)
        error('MESS:control_data', 'eqn.B is not defined or corrupted');
    end
    m = size(eqn.B, 2);
end

% make sure the first right hand side is dense so that the resulting factor
% is densly stored.
if isfield(eqn,'G')&&issparse(eqn.G), eqn.G=full(eqn.G); end
if isfield(eqn,'B')&&issparse(eqn.B), eqn.B=full(eqn.B); end
if isfield(eqn,'C')&&issparse(eqn.C), eqn.C=full(eqn.C); end
if isfield(eqn,'U')&&issparse(eqn.U), eqn.U=full(eqn.U); end
if isfield(eqn,'V')&&issparse(eqn.V), eqn.V=full(eqn.V); end

% check whether LDL^T formulation should be used
if ~isfield(opts.adi, 'LDL_T'), opts.adi.LDL_T = 0; end
% check for or set proper right hand side in eqn.G
if opts.adi.LDL_T
    % RHS has form G * S * G'
    % Solution has form L * D * L' with D Kronecker product of out.D and S
    % D is not computed explicitly
    % L is stored in Z if computed (opts.adi.computeZ)
    % G (eqn.G) and S (eqn.S) need to be given
    if ~isfield(eqn, 'G') || ~isnumeric(eqn.G)
        error('MESS:control_data', 'eqn.G is not defined or corrupted');
    end
    if ~isfield(eqn, 'S') || ~isnumeric(eqn.S)
        error('MESS:control_data', 'eqn.S is not defined or corrupted');
    end
    % init solution factor D
    out.D = zeros(opts.adi.maxiter, 1);
else
    if ~isfield(eqn, 'G')
        if eqn.type == 'N'
            eqn.G = eqn.B;
        else
            eqn.G = eqn.C';
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for feedback and shift matrices appearing inside
% Newton, BDF and Rosenbrock type methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts,'rosenbrock'), opts.rosenbrock=[]; end
if isstruct(opts.rosenbrock)&&isfield(opts.rosenbrock,'tau')
    rosenbrock = 1;
else
    rosenbrock = 0;
end
if ~isfield(opts,'bdf'), opts.bdf=[]; end
if isstruct(opts.bdf) && isfield(opts.bdf, 'tau') && isfield(opts.bdf, 'beta')
    bdf = 1;
else
    bdf = 0;
end

if ~isfield(eqn,'U') || isempty(eqn.U) || ~isfield(eqn,'V') || isempty(eqn.V)
    eqn.haveUV=0;
else
    if isnumeric(eqn.U) && isnumeric(eqn.V) && ...
            size(eqn.U,1)==size(eqn.V,1) && size(eqn.U,2)==size(eqn.V,2)
        eqn.haveUV=1;
    else
        error('MESS:SMW','Inappropriate data of low rank updated opertor (eqn.U and eqn.V)');
    end
end
% Get relevant sizes of right hand side and shift vector.
k=size(eqn.G,2);
l=length(opts.adi.shifts.p);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for projection triggers and residual control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project=isfield(opts.adi,'projection') && isstruct(opts.adi.projection) &&...
    isfield(opts.adi.projection,'freq')&&isnumeric(opts.adi.projection.freq)...
    &&opts.adi.projection.freq;

if project && (~isfield(opts.adi,'res') || ~isstruct(opts.adi.res))
    warning('MESS:control_data',...
        'Found empty residual control parameters. Falling back to defaults.')
    opts.adi.res=[];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether we want to compute Z or rather accumulate K in the Newton
% context for AREs, or both, e.g., in inexact Newton contexts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accumulation of K or DeltaK is helpful in inexact Newton and implicit
% Newton settings
if isfield(opts.adi,'accumulateK')&&opts.adi.accumulateK
    if eqn.type == 'T'
        out.Knew = zeros(size(eqn.B));
    else
        out.Knew = zeros([size(eqn.C, 2), size(eqn.C, 1)]);
    end
else
    opts.adi.accumulateK=0;
end
if isfield(opts.adi,'accumulateDeltaK')&&opts.adi.accumulateDeltaK
    if isfield(eqn, 'V') && ~isempty(eqn.V)  ...
            && ~(isfield(eqn, 'setUV') && eqn.setUV)
            out.DeltaK = -eqn.V; % eqn.setUV is false and in eqn.V is only K
    elseif isfield(eqn, 'V') && ~isempty(eqn.V) && (size(eqn.V, 2) > m)
            out.DeltaK = -eqn.V( : , end - m + 1 : end); % eqn.setUV is 
            % true and K is only in the last columns of eqn.V
    else % eqn.V is empty or in eqn.V is only the initial V from eqn.setUV
        % and K = []
        if eqn.type == 'T'
            out.DeltaK = zeros(size(eqn.B));
        else
            out.DeltaK = zeros([size(eqn.C, 2), size(eqn.C, 1)]);
        end
    end
else
    opts.adi.accumulateDeltaK=0;
end
if ~isfield(opts.adi, 'computeZ') || project
    opts.adi.computeZ = 1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All checks done. Here comes the real work!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize required usf for multiplication with E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eqn.haveE, [eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper); end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = [];
if opts.adi.restol, res=zeros(1, opts.adi.maxiter); else res=[]; end
if opts.adi.rctol, rc=zeros(1, opts.adi.maxiter); else rc=[]; end

i = 1;
i_shift = 1;

[eqn, opts, oper] = oper.init_res_pre(eqn, opts, oper);
[W, res0] = oper.init_res(eqn, opts, eqn.G);

if opts.adi.shifts.info
    fprintf('ADI Shifts:\n');
    opts.adi.shifts.p
end
out.linesearch = 0;
out.restart = 0;
if isfield(opts, 'nm') && isfield(opts.nm, 'accumulateRes') && ...
        opts.nm.accumulateRes && isfield(opts.nm, 'res0')
    outer_res = zeros(1, opts.adi.maxiter);
    res0 = opts.nm.res0;
else
    outer_res = [];
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while i<opts.adi.maxiter+1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check whether shifts need to be updated
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i_shift > l
        i_shift = 1;
        if opts.adi.LDL_T
            [ opts, l ] = mess_get_projection_shifts( eqn, opts, oper, ...
                Z, W, out.D(1 : i - 1));
        else
            [ opts, l ] = mess_get_projection_shifts( eqn, opts, oper, Z, W);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get current shift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pc = opts.adi.shifts.p(i_shift);
    projected=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform the actual step computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if bdf
        [ V, eqn, opts, oper ] = mess_solve_shifted_system_BDF(eqn, opts, oper, pc, W);
    elseif rosenbrock
        [ V, eqn, opts, oper ] = mess_solve_shifted_system_Rosenbrock(eqn, opts, oper, pc, W);
    else
        [ V, eqn, opts, oper ] = mess_solve_shifted_system(eqn, opts, oper, pc, W);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update low rank solution factor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isreal(pc)
        % just update the factor
        V=real(V);
        if opts.adi.computeZ
            if opts.adi.LDL_T
                Z = [Z, V];
                out.D(i) = -2 * pc;
            else
                Z = [Z, sqrt(-2*pc)*V];
            end
        end
        % update low rank residual
        if eqn.haveE
            EV = oper.mul_E(eqn, opts,eqn.type,V,'N');
            W = W - 2 * pc * EV;
        else
            W = W - 2 * pc * V;
        end
        if opts.adi.accumulateK
            if eqn.haveE
                if eqn.type == 'T'
                    if opts.adi.LDL_T
                        out.Knew=out.Knew+(oper.mul_E(eqn, opts,eqn.type,V,'N'))*...
                            ((2*(-pc) * eqn.S)*(V'*eqn.B));
                    else
                        out.Knew=out.Knew+(oper.mul_E(eqn, opts,eqn.type,V,'N'))*...
                            ((2*(-pc))*(V'*eqn.B));
                    end
                else
                    if opts.adi.LDL_T
                        out.Knew=out.Knew+(oper.mul_E(eqn, opts,eqn.type,V,'N'))*...
                            ((2*(-pc) * eqn.S)*(eqn.C * V)');
                    else
                        out.Knew=out.Knew+(oper.mul_E(eqn, opts,eqn.type,V,'N'))*...
                            ((2*(-pc))*(eqn.C * V)');
                    end
                end
            else
                if eqn.type == 'T'
                    if opts.adi.LDL_T
                        out.Knew=out.Knew+V*((2*(-pc) * eqn.S)*(V'*eqn.B));
                    else
                        out.Knew=out.Knew+V*((2*(-pc))*(V'*eqn.B));
                    end
                else
                    if opts.adi.LDL_T
                        out.Knew=out.Knew+V*((2*(-pc) * eqn.S)*(eqn.C * V)');
                    else
                        out.Knew=out.Knew+V*((2*(-pc))*(eqn.C * V)');
                    end
                end
            end
        end
        if opts.adi.accumulateDeltaK
            if eqn.haveE
                if eqn.type == 'T'
                    if opts.adi.LDL_T
                        out.DeltaK = out.DeltaK+(oper.mul_E(eqn, opts,eqn.type,V,'N'))*...
                            ((2*(-pc) * eqn.S)*(V'*eqn.B));
                    else
                        out.DeltaK = out.DeltaK+(oper.mul_E(eqn, opts,eqn.type,V,'N'))*...
                            ((2*(-pc))*(V'*eqn.B));
                    end
                else
                    if opts.adi.LDL_T
                        out.DeltaK = out.DeltaK+(oper.mul_E(eqn, opts,eqn.type,V,'N'))*...
                            ((2*(-pc) * eqn.S)*(eqn.C * V)');
                    else
                        out.DeltaK = out.DeltaK+(oper.mul_E(eqn, opts,eqn.type,V,'N'))*...
                            ((2*(-pc))*(eqn.C * V)');
                    end
                end
            else
                if eqn.type == 'T'
                    if opts.adi.LDL_T
                        out.DeltaK = out.DeltaK+V*((2*(-pc) * eqn.S)*(V'*eqn.B));
                    else
                        out.DeltaK = out.DeltaK+V*((2*(-pc))*(V'*eqn.B));
                    end
                else
                    if opts.adi.LDL_T
                        out.DeltaK = out.DeltaK+V*((2*(-pc) * eqn.S)*(eqn.C * V)');
                    else
                        out.DeltaK = out.DeltaK+V*((2*(-pc))*(eqn.C * V)');
                    end
                end
            end
        end
    else
        % perform a double step with the known solution for the conjugate shift
        a=2*sqrt(-real(pc)); b=real(pc)/imag(pc);
        V1=a*(real(V)+b*imag(V));
        V2=(a*sqrt(b*b+1))*imag(V);
        if opts.adi.computeZ
            if opts.adi.LDL_T
                Z = [Z, (sqrt(2) / a) * V1, (sqrt(2) / a) * V2];
                out.D(i : i + 1) = -2 * real(pc);
            else
                Z = [Z, V1, V2];
            end
        end
        if opts.adi.accumulateK
            if eqn.haveE
                if eqn.type == 'T'
                    if opts.adi.LDL_T
                        out.Knew=out.Knew+oper.mul_E(eqn, opts,eqn.type,V1,'N')*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (V1'*eqn.B))...
                            +(oper.mul_E(eqn, opts,eqn.type,V2,'N'))*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (V2'*eqn.B));
                    else
                        out.Knew=out.Knew+oper.mul_E(eqn, opts,eqn.type,V1,'N')*(V1'*eqn.B)...
                            +(oper.mul_E(eqn, opts,eqn.type,V2,'N'))*(V2'*eqn.B);
                    end
                else
                    if opts.adi.LDL_T
                        out.Knew=out.Knew+oper.mul_E(eqn, opts,eqn.type,V1,'N')*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (eqn.C * V1)')...
                            +(oper.mul_E(eqn, opts,eqn.type,V2,'N'))*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (eqn.C * V2)');
                    else
                        out.Knew=out.Knew+oper.mul_E(eqn, opts,eqn.type,V1,'N')*(eqn.C * V1)'...
                            +(oper.mul_E(eqn, opts,eqn.type,V2,'N'))*(eqn.C * V2)';
                    end
                end
            else
                if eqn.type == 'T'
                    if opts.adi.LDL_T
                        out.Knew=out.Knew+V1*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (V1'*eqn.B))...
                            +V2*((2 / (a * a) * -2 * real(pc)) * eqn.S * (V2'*eqn.B));
                    else
                        out.Knew=out.Knew+V1*(V1'*eqn.B)+V2*(V2'*eqn.B);
                    end
                else
                    if opts.adi.LDL_T
                        out.Knew=out.Knew+V1*...
                             ((2 / (a * a) * -2 * real(pc)) * eqn.S * (eqn.C * V1)')...
                             +V2* ((2 / (a * a) * -2 * real(pc)) * eqn.S * (eqn.C * V2)');
                    else
                        out.Knew=out.Knew+V1*(eqn.C * V1)'+V2*(eqn.C * V2)';
                    end
                end
                
            end
        end
        if opts.adi.accumulateDeltaK
            if eqn.haveE
                if eqn.type == 'T'
                    if opts.adi.LDL_T
                        out.DeltaK=out.DeltaK+oper.mul_E(eqn, opts,eqn.type,V1,'N')*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (V1'*eqn.B))...
                            +(oper.mul_E(eqn, opts,eqn.type,V2,'N'))*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (V2'*eqn.B));
                    else
                        out.DeltaK=out.DeltaK+oper.mul_E(eqn, opts,eqn.type,V1,'N')*(V1'*eqn.B)...
                            +(oper.mul_E(eqn, opts,eqn.type,V2,'N'))*(V2'*eqn.B);
                    end
                else
                    if opts.adi.LDL_T
                        out.DeltaK=out.DeltaK+oper.mul_E(eqn, opts,eqn.type,V1,'N')*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (eqn.C * V1)')...
                            +(oper.mul_E(eqn, opts,eqn.type,V2,'N'))*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (eqn.C * V2)');
                    else
                        out.DeltaK=out.DeltaK+oper.mul_E(eqn, opts,eqn.type,V1,'N')*(eqn.C * V1)'...
                            +(oper.mul_E(eqn, opts,eqn.type,V2,'N'))*(eqn.C * V2)';
                    end
                end
            else
                if eqn.type == 'T'
                    if opts.adi.LDL_T
                        out.DeltaK=out.DeltaK+V1*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (V1'*eqn.B))...
                            +V2*((2 / (a * a) * -2 * real(pc)) * eqn.S * (V2'*eqn.B));
                    else
                        out.DeltaK=out.DeltaK+V1*(V1'*eqn.B)+V2*(V2'*eqn.B);
                    end
                else
                    if opts.adi.LDL_T
                        out.DeltaK=out.DeltaK+V1*...
                            ((2 / (a * a) * -2 * real(pc)) * eqn.S * (eqn.C * V1)')...
                            +V2*((2 / (a * a) * -2 * real(pc)) * eqn.S * (eqn.C * V2)');
                    else
                        out.DeltaK=out.DeltaK+V1*(eqn.C * V1)'+V2*(eqn.C * V2)';
                    end
                end
                
            end
        end
        % update low rank residual for double step
        if eqn.haveE
            EV = oper.mul_E(eqn, opts,eqn.type,real(V)+b*imag(V),'N');
            W = W - 4 * real(pc) * EV;
        else
            W = W + a * V1;
        end
        i=i+1;
        i_shift = i_shift + 1;
        if ~isempty(outer_res)
            if i > 2
                outer_res(i - 1) = outer_res(i - 2);
            else
                outer_res(i - 1) = opts.nm.res0;
            end
        end
        if i > 2 && opts.adi.restol~=0
            res(i - 1) = res(i - 2);
        else
            res(i - 1) = res0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform projection update if desired
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if project && (~mod(i, opts.adi.projection.freq))
        projected=1;
        Z=mess_galerkin_projection_acceleration(Z,'lyapunov',eqn,oper,opts);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.adi.restol
        if ~projected
            if opts.adi.LDL_T
                res(i) = max(abs(eig(W' * W * eqn.S))) / res0;
            else
                res(i) = norm(W' * W, opts.adi.norm) / res0;
            end
            if ~isempty(outer_res)
                outer_res(i) = riccati_LR(W, out.DeltaK, opts, eqn.S, []) ...
                    / opts.nm.res0;
            end
        else
            res(i)=res2_norms(Z,'lyapunov',eqn,opts,oper)/res0;
            if ~isempty(outer_res)
                outer_res(i) = res2_norms(Z, 'riccati', eqn, opts, oper) ...
                    / opts.nm.res0;
            end
        end
    end
    if opts.adi.rctol
        if ~projected
            if isreal(pc)
                nrmV=-2 * pc * sum(sum(V.^2));
            else % complex double step means 2 blocks added
                nrmV=sum(sum([V1, V2].^2));
            end
            nrmZ=nrmZ+nrmV;
        else
            if isreal(pc)
                nrmV=sum(sum(Z(:,(end-k+1):end).^2));
            else % complex double step means 2 blocks added
                nrmV=sum(sum(Z(:,(end-2*k+1):end).^2));
            end
            nrmZ=sum(sum(Z.^2));
        end
        rc(i) = sqrt(nrmV/nrmZ);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print status information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.adi.info
        if opts.adi.rctol&&opts.adi.restol
            fprintf(1,['ADI step: %4d normalized residual: %e relative change ' ...
                'in Z: %e\n'],i,res(i),rc(i));
        elseif opts.adi.restol
            fprintf(1,'ADI step: %4d normalized residual: %e \n',i,res(i));
        elseif opts.adi.rctol
            fprintf(1,['ADI step: %4d relative change ' ...
                'in Z: %e\n'],i,rc(i));
        end
        if ~isempty(outer_res)
            fprintf(1, '\t\tnormalized outer residual: %e\n', outer_res(i));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ opts, out, stop ] = prepare_next_adi_iteration( opts, out, res, ...
        rc, outer_res, i);
    if stop
        break
    end
    i=i+1;
    i_shift = i_shift + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print outer tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.adi.info && opts.adi.inexact
    fprintf(1, '\nouter tolerance: %e\n', opts.adi.outer_tol);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare output arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.niter = i - (i > opts.adi.maxiter);
if opts.adi.restol, out.res=res(1:out.niter); end
if opts.adi.rctol, out.rc=rc(1:out.niter); end
if opts.adi.LDL_T, out.D = out.D(1 : out.niter); end
out.res_fact = W;
if ~opts.adi.computeZ
    Z = [];
end
if ~isempty(outer_res)
    out.Riccati_res = outer_res(out.niter);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalize required usf for multiplication with E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eqn.haveE, [eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper); end
[eqn, opts, oper] = oper.init_res_post(eqn, opts, oper);
