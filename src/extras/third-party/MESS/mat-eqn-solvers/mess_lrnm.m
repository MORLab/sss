function [Z, out]=mess_lrnm(eqn, opts, oper)
%% function [Z,out]=mess_lrnm(eqn,opts, oper)
%
% Solve continuous-time Riccati equations with sparse coefficients with
% Newton's method (NM)
%   eqn.type = 'N' -> A*Z*Z'*E' + E*Z*Z'*A' - E*Z*Z'*C'*C*Z*Z'*E' + B*B' = 0 (N)
%   eqn.type = 'T' -> A'*Z*Z'*E + E'*Z*Z'*A - E'*Z*Z'*B*B'*Z*Z'*E + C'*C = 0 (T)
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
% Input fields in struct eqn:
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
%   eqn.K0      dense (n x m1) matrix, initial K
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
%   Depending on the operator chosen by the operatormanager, additional
%   fields may be needed. For the "default", e.g., eqn.A_ and eqn.E_ hold
%   the A and E matrices. For the second order types these are given
%   implicitly by the M, D, K matrices stored in eqn.M_, eqn.D_ and eqn.K_,
%   respectively.
%
% Input fields in struct opts:
%   opts.nm.maxiter             possible  values: integer > 0
%                               maximum NM iteration number
%                               (optional)
%
%   opts.nm.restol              possible  values: scalar >= 0
%                               stopping tolerance for the relative NM
%                               residual norm; if restol = 0 the relative
%                               residual norm is not evaluated
%                               (optional)
%
%   opts.nm.rctol               possible  values: scalar >= 0
%                               stopping tolerance for the relative
%                               change of the NM solutiuon Z; 
%                               if restol = 0 the relative
%                               change is not evaluated
%                               (optional)
%
%   opts.nm.info                possible  values: 0, 1, false, true
%                               turn on (1) or off (0) the status output in
%                               every NM iteration step
%                               (optional)
%
%   opts.nm.norm                possible  values: 2, 'fro'
%                               use 2-norm (2) or Frobenius norm ('fro') to
%                               compute residual and relative change norms
%                               (optional)
%
%   opts.nm.accumulateRes       possible  values: 0, 1, false, true
%                               accumulate the relative NM residual norm
%                               during the inner ADI iteration
%                               (optional)
%
%   opts.nm.linesearch          possible  values: 0, 1, false, true
%                               if tuned of (0) NM makes full steps; if
%                               turned on (1) a step size 0<=lambda<=2 is
%                               computed
%                               (optional)
%
%   opts.nm.inexact             possible  values: 0, false, 'linear', 
%                               'superlinear', 'quadratic'
%                               the inner ADI uses an adaptive relative ADI
%                               residual norm; with 
%                               ||R||: relative NM residual norm,
%                               tau: opts.nm.tau,
%                               j: NM iteration intex;
%                               'linear': tau * ||R||,  
%                               'superlinear':  ||R|| / (j^3 + 1),
%                               'quadratic': tau / sqrt(||R||)
%                               (optional)
%                               
%   opts.nm.tau                 possible  values: scalar >= 0
%                               factor for inexact inner ADI iteration
%                               tolerance
%                               (optional)
%
%   opts.nm.LDL_T               possible  values: 0, 1, false, true
%                               solution is represented as L*D*L' instead
%                               of Z*Z'; eqn.S needs to be given
%                               (optional)
% 
%   opts.nm.projection.freq     possible  values: integer > 0
%                               frequency of the usage of galerkin
%                               projection acceleration in NM
%                               (optional)
%
%   opts.nm.projection.ortho    possible  values: 0, 1, false, true
%                               implicit (0) or explicit (1) 
%                               orthogonalization of Z in NM
%                               (optional)
%
%   opts.nm.projection.meth     possible  values: 'lyapchol', 
%                               'lyap_sgn_fac', 'lyap', 'lyapunov',
%                               'lyap2solve'
%                               projection method in NM
%                               (optional)
%
%   opts.adi.maxiter            possible  values: integer > 0
%                               maximum ADI iteration number
%                               (optional)
%
%   opts.adi.restol             possible  values: scalar >= 0
%                               stopping tolerance for the relative
%                               ADI residual norm; if restol = 0 the 
%                               relative residual norm is not evaluated
%                               (optional)
%
%   opts.adi.rctol              possible  values: scalar >= 0
%                               stopping tolerance for the relative
%                               change of the ADI solutiuon Z; 
%                               if restol = 0 the relative
%                               change is not evaluated
%                               (optional)
%
%   opts.adi.info               possible  values: 0, 1, false, true
%                               turn on (1) or off (0) the status output in
%                               every ADI iteration step
%                               (optional)
%
%   opts.adi.norm               possible  values: 2, 'fro'
%                               use 2-norm (2) or Frobenius norm ('fro') to
%                               compute residual and relative change norms;
%                               must be the same as opts.nm.norm
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
%                               ADI iteration
%                               (optional)
%
%   opts.adi.accumulateDeltaK   possible  values: 0, 1, false, true
%                               accumulate the update DeltaK of the 
%                               feedback matrix K during the ADI iteration
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
%                               projection acceleration in ADI
%                               (optional)
%
%   opts.adi.projection.ortho   possible  values: 0, 1, false, true
%                               implicit (0) or explicit (1) 
%                               orthogonalization of Z in ADI
%                               (optional)
%
%   opts.adi.projection.meth    possible  values: 'lyapchol', 
%                               'lyap_sgn_fac', 'lyap', 'lyapunov',
%                               'lyap2solve'
%                               projection method in ADI
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
% For LDL^T formulation use opts.nm.LDL_T = 1:
%     RHS of Lyapunov Eq. has form G * S * G'
%     Solution Lyapunov Eq. has form L * D * L' 
%     with D Kronecker product of adiout.D and S
%     D is not computed explicitly
%     L is stored in Z if computed (opts.adi.computeZ)
%     S (eqn.S) needs to be given
%
% Output fields in struct out:
%   out.adi             struct with the output of the last ADI iteration
%
%   out.nm.niter        number of NM iterations
%
%   out.nm.K            feedback matrix K
%
%   out.nm.D            scalar factor vector for LDL^T formulation
%                       use kron(diag(out.nm.D),eqn.S) to build D
%                       (opts.adi.LDL_T = 1)
%
%   out.nm.S            matrix S for LDL^T formulation
%                       use kron(diag(out.nm.D),out.nm.S) to build D
%                       (opts.adi.LDL_T = 1)
%
%   out.nm.res          array of relative NM residual norms
%
%   out.nm.rc           array of relative NM change norms
%
%   uses operatorfunctions init and mul_E, mul_E_pre, mul_E_post directly
%   and indirectly:
%     size, init, init_res, sol_ApE and mul_E in mess_lradi;
%     size in mess_para
%     size, sol_A, mul_A, sol_E, mul_E in mess_arn
%     mul_A, mul_E in mess_galerkin_projection_acceleration
%     mul_A, mul_E in ricatti
%

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for ADI Control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(opts,'adi') || ~isstruct(opts.adi)
    error('MESS:control_data','ADI control structure opts.ADI missing.');
end % Single fields are checked below or inside mess_lradi
if ~isfield(opts.adi,'computeZ'),  opts.adi.computeZ=1; end
if ~isfield(opts.adi,'accumulateK'),  opts.adi.accumulateK=0; end
if ~isfield(opts.adi,'accumulateDeltaK'),  opts.adi.accumulateDeltaK=0; end
if ~(opts.adi.computeZ || opts.adi.accumulateK || opts.adi.accumulateDeltaK)
    warning('MESS:control_data', ...
        ['Either opts.adi.accumulateK or opts.adi.computeZ or ', ...
        'opts.adi.accumulateDeltaK must be 1.']);
    opts.adi.accumulateDeltaK = 1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for shift parameter structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts.adi,'shifts') || ~isstruct(opts.adi.shifts)
    error('MESS:control_data','shift parameter control structure missing.');
else
    if~isfield(opts.adi.shifts,'l0')||...
            ~isfield(opts.adi.shifts,'kp')||...
            ~isfield(opts.adi.shifts,'km')
        error('MESS:shifts', ...
            'Incomplete input parameters for shift computation.');
    end
end
if ~isfield(opts.adi.shifts,'period')
    opts.adi.shifts.period=1;
end
if ~isfield(opts.adi.shifts,'method')
    opts.adi.shifts.method='heur';
end
if strcmp(opts.adi.shifts.method, 'wachspress')
    if ~isfield(opts.adi.shifts,'wachspress')
        opts.adi.shifts.wachspress='T';
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for Newton control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts,'nm') || ~isstruct(opts.nm)
    error('MESS:control_data','Newton control structure opts.nm missing.');
else
    if ~isfield(opts.nm,'maxiter')||~isnumeric(opts.nm.maxiter)
        warning('MESS:control_data',...
            'Missing or Corrupted maxiter field. Switching to default.');
        opts.nm.maxiter=20;
    end
    
    if ~isfield(opts.nm,'rctol')||~isnumeric(opts.nm.rctol)
        warning('MESS:control_data',...
            'Missing or Corrupted rctol field. Switching to default.');
        opts.nm.rctol=0;
    end
    
    if ~isfield(opts.nm,'restol')||~isnumeric(opts.nm.restol)
        warning('MESS:control_data',...
            'Missing or Corrupted restol field. Switching to default.');
        opts.nm.restol=0;
    end
    
    if ~isfield(opts.nm,'accumulateRes')||~isnumeric(opts.nm.accumulateRes)
        warning('MESS:control_data', ...
            'Missing or Corrupted accumulateRes field. Switching to default.');
        opts.nm.accumulateRes = 0;
    end
    if opts.nm.accumulateRes
        % need DeltaK
        opts.adi.accumulateDeltaK = 1;
    end
    if ~isfield(opts.nm, 'linesearch') || ~isnumeric(opts.nm.linesearch)
        warning('MESS:control_data', ...
            'Missing or Corrupted linesearch field. Switching to default.');
        opts.nm.linesearch = 0;
    end
    if opts.nm.linesearch
        % need DeltaK
        opts.adi.accumulateDeltaK = 1;
        % need restol
        if ~opts.nm.restol
            opts.nm.restol = 1e-16;
        end
        alpha = 1e-4;
    end
    if ~isfield(opts.nm,'inexact'),  opts.nm.inexact=0; end
    if opts.nm.inexact
        if ~isfield(opts.nm,'tau'),  opts.nm.tau=1; end
        if ~isfield(opts.adi,'restol'),  opts.adi.restol=1e-16; end
        opts.nm.accumulateRes = 1;
        opts.adi.inexact = 1;
        opts.adi.accumulateDeltaK = 1;
    else
        opts.adi.inexact = 0;
    end
    if ~isfield(opts.nm, 'norm') || ...
            (~strcmp(opts.nm.norm, 'fro') && ...
            (~isnumeric(opts.nm.norm) ||  opts.nm.norm ~= 2))
        warning('MESS:control_data', ...
            'Missing or Corrupted norm field. Switching to default.');
        opts.nm.norm = 'fro';
        opts.adi.norm = 'fro';
    else
        opts.adi.norm = opts.nm.norm;
    end
    
    if ~isfield(opts.nm,'info')||~isnumeric(opts.nm.restol)
        warning('MESS:control_data',...
            'Missing or Corrupted info field. Switching to default.');
        opts.nm.info=0;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
[eqn, erg] = oper.init(eqn, opts, 'A','E');
if ~erg
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end

if ~isfield(eqn, 'B') || ~isnumeric(eqn.B)
    error('MESS:control_data', 'eqn.B is not defined or corrupted');
end

if ~isfield(eqn, 'C') || ~isnumeric(eqn.C)
    error('MESS:control_data', 'eqn.C is not defined or corrupted');
end

% make sure the first right hand side is dense so that the resulting factor
% is densly stored.
if issparse(eqn.C), eqn.C = full(eqn.C); end
if issparse(eqn.B), eqn.B = full(eqn.B); end

if ~isfield(eqn, 'type')
    eqn.type = 'N';
    warning('MESS:control_data',['Unable to determine type of equation.'...
        'Falling back to type ''N''']);
elseif (eqn.type ~= 'N') && (eqn.type ~= 'T')
    error('MESS:equation_type', 'Equation type must be either ''T'' or ''N''');
end

if eqn.type == 'T'
    p=size(eqn.C,1); %number of outputs
    m=size(eqn.B,2); %number of inputs
else
    m=size(eqn.C,1); %number of outputs
    p=size(eqn.B,2); %number of inputs
end
% check whether LDL^T formulation should be used
if ~isfield(opts.nm, 'LDL_T'), opts.nm.LDL_T = 0; end
% check for or set proper right hand side
if opts.nm.LDL_T
    % RHS of Lyapunov Eq. has form G * S * G'
    % Solution Lyapunov Eq. has form L * D * L' 
    % with D Kronecker product of adiout.D and S
    % D is not computed explicitly
    % L is stored in Z if computed (opts.adi.computeZ)
    % S (eqn.S) need to be given
    if ~isfield(eqn, 'S') || ~isnumeric(eqn.S)
        error('MESS:control_data', 'eqn.S is not defined or corrupted');
    end
    if isfield(opts, 'bdf') && isstruct(opts.bdf)
        if ~isfield(opts.bdf, 'tau') || ~isnumeric(opts.bdf.tau)
            error('MESS:control_data', 'opts.bdf.tau is not defined or corrupted');
        end
        if ~isfield(opts.bdf, 'beta') || ~isnumeric(opts.bdf.beta)
            error('MESS:control_data', 'opts.bdf.beta is not defined or corrupted');
        end
        opts.adi.LDL_T = 1;
        tau_beta = opts.bdf.tau * opts.bdf.beta;
        bdf = 1;
    else
        bdf = 0;
        opts.adi.LDL_T = 1;
        tau_beta = 1;
    end
else
    bdf = 0;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% System data for Riccati iteration                                   % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ~isfield( eqn, 'setUV' )                                             
    eqn.setUV = 0;  
    s = 0;
else                                                                    
    if eqn.setUV  
        if opts.nm.LDL_T || (isfield(opts.adi, 'LDL_T') && opts.adi.LDL_T)
            error('MESS:control_data', ...
                'LDL_T formulation is not compatible with eqn.setUV option');
        end
        if issparse( eqn.V )                                            
            eqn.V = full( eqn.V );                                      
        end                                                             
        if issparse( eqn.U )                                            
            eqn.U = full( eqn.U );                                      
        end                                                             
        if ~isfield( eqn, 'UVsize' )                                    
            eqn.UVsize = size( eqn.U, 2 );                              
        end
        s = eqn.UVsize;   
    else
        s = 0;
    end                                                                 
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if we have an initial stabilizing feedback and initialize storage
% for the computed feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(eqn,'K0')
    if eqn.setUV
        eqn.V = [ eqn.V, eqn.K0 ];
    else
        eqn.V = eqn.K0;
    end
    if eqn.type == 'T'
        eqn.U( : , s + 1 : s + m) = eqn.B;
    else
        eqn.U( : , s + 1 : s + m) = eqn.C';
    end
    eqn.haveUV=1;
else
    if eqn.setUV
        eqn.haveUV = 1;
    else
        eqn.V      = [];
        eqn.haveUV = 0;
    end
end
if eqn.type == 'T';
    if bdf
        eqn.U = opts.bdf.tau * opts.bdf.beta * eqn.B;
    end
else
    if bdf
        eqn.U = opts.bdf.tau * opts.bdf.beta * eqn.C';
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for projection triggers and residual control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project=isfield(opts.nm,'projection')&&...
    isfield(opts.nm.projection,'freq')&&isnumeric(opts.nm.projection.freq)...
    &&opts.nm.projection.freq;
if project
    opts.adi.computeZ = 1;
    opts.nm.norm = 2;
    opts.adi.norm = 2;
end

% we need to use iterative residual computation. Let's add the
% corresponding control structure in case it does not already exist.
if project && ~isfield(opts.nm,'res')
    warning('MESS:control_data',...
        'Found empty residual control parameters. Falling back to defaults.')
    opts.nm.res=[];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All checks done. Here comes the real work!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if changes are made here these have to be repeated at the restart part
% from lineserach further below!!!
if opts.nm.restol, res=zeros(opts.nm.maxiter,1); else res=[]; end
if opts.nm.rctol, rc=zeros(opts.nm.maxiter,1); else rc=[]; end

if eqn.type == 'T'
    eqn.G = eqn.C';
else
    eqn.G = eqn.B;
end

if eqn.setUV
    eqn.G = [ eqn.G, eqn.V( : , s + 1 : end) ];
else
    eqn.G = [ eqn.G, eqn.V ];
end
if opts.nm.LDL_T
    res0 = riccati_LR(eqn.G, [], opts, eqn.S, []);
    S = eqn.S;
    if bdf
        eqn.S = blkdiag(eqn.S, tau_beta * eye(size(eqn.V, 2)));
    end
else
    res0=norm(eqn.G'*eqn.G, opts.nm.norm);
    eqn.S = [];
end
if opts.nm.accumulateRes
    opts.nm.res0 = res0;
end
if opts.adi.inexact
    switch opts.nm.inexact
        case 'linear'
            opts.adi.outer_tol = opts.nm.tau * res0;
        case 'superlinear'
            opts.adi.outer_tol = 0.5 * res0;
        case 'quadratic'
            if res0 > 1
                opts.adi.outer_tol = opts.nm.tau / sqrt(res0);
            else
                opts.adi.outer_tol = opts.nm.tau * res0 * res0;
            end
        otherwise
            error('MESS:inexact', ...
                ['inexact must be 0, ''linear'', ''superlinear''', ...
                ' or ''quadratic''']);
    end
    opts.adi.outer_tol = max(opts.adi.outer_tol, opts.adi.restol);
end
if opts.nm.linesearch
    linesearch = 0;
    W_old = eqn.G;
    if ~isfield(eqn,'K0')
        DeltaK_old = [];
    else
        %         DeltaK_old = -eqn.V( : , s + 1 : end); % careful with setUV
        %         option!
        DeltaK_old = [];
    end
    if opts.nm.LDL_T
        S_old = S;
    end
end
restarted = 0;
already_restarted = 0;
% if changes are made here these have to be repeated at the restart part
% from lineserach further below!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eqn,opts,oper]=oper.mul_E_pre(eqn,opts,oper);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=1;
while j<=opts.nm.maxiter
    projected=0;
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute new ADI shifts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~mod(j - 1,opts.adi.shifts.period)
        opts.adi.shifts.p=mess_para(eqn,opts,oper);
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform the actual step computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % form right hand side factor
    if size(eqn.V, 2) > s % only possible in first step with no initial
        % stabilizing feedback
        eqn.G( : , p + 1 : p + m) = eqn.V( : , end - m + 1 : end);
        if opts.nm.LDL_T
            eqn.S = blkdiag(S, tau_beta * eye(size(eqn.V, 2)));
        end
    end
    % solve the Lyapunov equation
    [Z, adiout,eqn,opts,oper]=mess_lradi(eqn,opts, oper);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Restart Newton iteration with exact ADI iteration if necessary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if adiout.restart && opts.nm.inexact(1)
        if already_restarted
            error('MESS:lrnm', ...
                        'Newton iteration with line search failed to converge.');
        else
            %restart Newton iteration
            warning('MESS:lrnm', ...
                ['Newton iteration needs to be restartet because ', ...
                'of divergence. Continuing with exact ADI iteration.']);
            if opts.nm.restol, res=[res; zeros(opts.nm.maxiter,1)]; else res=[]; end
            if opts.nm.rctol, rc=[rc; zeros(opts.nm.maxiter,1)]; else rc=[]; end
            opts.nm.maxiter = opts.nm.maxiter + j;
            if isfield(eqn,'K0')
                if eqn.setUV
                    eqn.V( : , (end - m + 1) : end) = eqn.K0;
                else
                    eqn.V = eqn.K0;
                end
                eqn.haveUV=1;
            else
                if eqn.setUV
                    eqn.haveUV = 1;
                else
                    eqn.V      = [];
                    eqn.haveUV = 0;
                end
            end
            if eqn.type == 'T'
                eqn.G = eqn.C';
            else
                eqn.G = eqn.B;
            end
            
            if eqn.setUV
                eqn.G = [ eqn.G, eqn.V( : , s + 1 : end) ];
            else
                eqn.G = [ eqn.G, eqn.V ];
            end
            if opts.nm.LDL_T
                res0 = riccati_LR(eqn.G, [], opts, eqn.S, []);
                S = eqn.S;
                if bdf
                    eqn.S = blkdiag(eqn.S, tau_beta * eye(size(eqn.V, 2)));
                end
            else
                res0=norm(eqn.G'*eqn.G, opts.nm.norm);
                eqn.S = [];
            end
            if opts.nm.accumulateRes
                opts.nm.res0 = res0;
            end
            opts.nm.inexact = 0;
            opts.adi.inexact = 0;
            if opts.nm.linesearch
                linesearch = 0;
                W_old = eqn.G;
                if ~isfield(eqn,'K0')
                    DeltaK_old = [];
                else
                    %         DeltaK_old = -eqn.V(:,end-m+1:end); % careful with setUV
                    %         option!
                    DeltaK_old = [];
                end
                if opts.nm.LDL_T
                    S_old = S;
                end
            end
            restarted = 1;
            continue
        end
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform projection update if desired
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if project&&(~mod(j,opts.nm.projection.freq))
        projected=1;
        Z=mess_galerkin_projection_acceleration(Z,'riccati',eqn,oper,opts);
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update feedback approximation if it has not been accumulated during 
    % the inner iteration or after projection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~opts.adi.accumulateK || projected
        if ~opts.adi.accumulateDeltaK || projected
            if eqn.haveE
                if opts.nm.LDL_T
                    if eqn.type == 'T'
                        adiout.Knew = oper.mul_E(eqn, opts,eqn.type,Z,'N')*...
                            (kron(diag(adiout.D), eqn.S) * (Z'*eqn.B));
                    else
                        adiout.Knew = oper.mul_E(eqn, opts,eqn.type,Z,'N')*...
                            (kron(diag(adiout.D), eqn.S) * (eqn.C * Z)');
                    end
                else
                    if eqn.type == 'T'
                        adiout.Knew = oper.mul_E(eqn, opts,eqn.type,Z,'N')*(Z'*eqn.B);
                    else
                        adiout.Knew = oper.mul_E(eqn, opts,eqn.type,Z,'N')*(eqn.C * Z)';
                    end
                end
            else
                if opts.nm.LDL_T
                    if eqn.type == 'T'
                        adiout.Knew=Z*(kron(diag(adiout.D), eqn.S) * (Z'*eqn.B));
                    else
                        adiout.Knew=Z*(kron(diag(adiout.D), eqn.S) * (eqn.C * Z)');
                    end
                else
                    if eqn.type == 'T'
                        adiout.Knew=Z*(Z'*eqn.B);
                    else
                        adiout.Knew=Z*(eqn.C * Z)';
                    end
                end
            end
        else
            if size(eqn.V, 2) == s
                adiout.Knew = adiout.DeltaK;
            else
                adiout.Knew = eqn.V(:,end-m+1:end) + adiout.DeltaK;
            end
        end
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.nm.restol||opts.nm.rctol
        if size(eqn.V, 2) > s
            V1=adiout.Knew-eqn.V(:,end-m+1:end);
        else
            V1=adiout.Knew;
        end
    end
    if opts.nm.restol
        if projected
            res(j)=res2_norms(Z,'riccati',eqn,opts,oper)/res0;
        else
            if opts.nm.accumulateRes
                res(j) = adiout.Riccati_res;
            else
                res(j)=norm([V1 adiout.res_fact], opts.nm.norm)^2/res0;
            end
        end
    end
    if opts.nm.linesearch && ~projected
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check whether line search is necessary
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ((j == 1) && (res(j) > 1)) || ((j > 1) && (res(j) > res(j - 1))) ...
                || adiout.linesearch
            linesearch = 1;
            % Compute Lambda
            lambda = exact_line_search(W_old, DeltaK_old, adiout.res_fact, ...
                adiout.DeltaK);
            if opts.nm.info
                fprintf(['\n\t\tUsing line search (res: %4d)\n',...
                    '\t\tlambda: %e\n'], res(j), lambda);
            end
            % Update K
            if size(eqn.V, 2) == s
                adiout.Knew = lambda * adiout.DeltaK;
            else
                adiout.Knew = eqn.V(:,end-m+1:end) + lambda * adiout.DeltaK;
            end
            % Update DeltaK and W
            if lambda <= 1
                adiout.DeltaK = [sqrt(1 - lambda) * DeltaK_old, ...
                    lambda * adiout.DeltaK];
                adiout.res_fact = [sqrt(1 - lambda) * W_old, ...
                    sqrt(lambda) * adiout.res_fact];
                if opts.nm.LDL_T
                    S_linesearch = blkdiag(S_old, eqn.S);
                    S_K = [];
                end
            else
                W_new = adiout.res_fact;
                DeltaK_new = adiout.DeltaK;
                sqrt_lambda = sqrt(lambda);
                adiout.DeltaK = [W_old, sqrt_lambda * W_new, ...
                    sqrt_lambda * DeltaK_old];
                adiout.res_fact = [DeltaK_old, sqrt_lambda * W_old, ...
                    lambda * DeltaK_new];
                if opts.nm.LDL_T
                    S_linesearch = blkdiag(eye(size(DeltaK_old, 2)), S_old, ...
                        eye(size(DeltaK_new, 2)));
                    S_K = blkdiag(S_old, eqn.S, eye(size(DeltaK_old, 2)));
                end
            end
            if ~opts.nm.LDL_T
                S_linesearch = eqn.S;
                S_K = [];
            end
            % Compute residual norm after line search
            res(j) = riccati_LR(adiout.res_fact, adiout.DeltaK, opts, ...
                S_linesearch, S_K) / res0;
            if ~restarted && (((j == 1) && (res(j) >= (1 - lambda * alpha))) || ...
                    ((j > 1) && (res(j) >= (1 - lambda * alpha) * res(j - 1))))
                % No sufficient decrease
                warning('MESS:lrnm', ...
                    ['Newton iteration with line search has', ...
                    ' insufficient decrease. ']);
                if opts.adi.inexact
                    % switch to exact ADI iteration
                    warning('MESS:lrnm', ...
                    'Switching to exact ADI iteration.');
                    opts.adi.inexact = 0;
                else
                    % Newton iteration with line search failed to converge, stop
                    error('MESS:lrnm', ...
                        'Newton iteration with line search failed to converge.');
                end
            end
        end
        % Keep ADI LR residual factor and DeltaK for next Newton iteration
        W_old = adiout.res_fact;
        DeltaK_old = adiout.DeltaK;
        if opts.nm.LDL_T
            if linesearch
                S_old = S_linesearch;
                % be careful here, adiout.res_fact, adiout.DeltaK don't fit
                % eqn.S anymore, till now they are not needed anywhere
                % together after this point
            else
                S_old = eqn.S;
            end
        end
        linesearch = 0;
    end
    
    if opts.nm.rctol
        if opts.adi.accumulateDeltaK
            rc(j) = norm(adiout.DeltaK,  opts.nm.norm) ...
                / norm(adiout.Knew,  opts.nm.norm);
        else
            rc(j)=norm(V1, opts.nm.norm)/norm(adiout.Knew, opts.nm.norm);
        end
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print status information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.nm.info
        if opts.nm.rctol&&opts.nm.restol
            fprintf(1,['\nNM step: %4d  normalized residual: %e\n'...
                '               relative change in K: %e\n'...
                '               number of ADI steps: %d \n\n'],...
                j,res(j),rc(j),adiout.niter);
        elseif opts.nm.restol
            fprintf(1,['\nNM step: %4d  normalized residual: %e \n\n'...
                '                number of ADI steps: %d \n\n']...
                ,j,res(j),adiout.niter);
        elseif opts.nm.rctol
            fprintf(1,['\nNM step: %4d  relative change in K: %e\n\n'...
                '               number of ADI steps: %d \n\n']...
                ,j,rc(j));
        end
    end
    if eqn.setUV
        eqn.V( : , s + 1 : s + m) = adiout.Knew;
    else
        eqn.V=adiout.Knew;
    end
    if ~isfield(eqn, 'U') || (size(eqn.U, 2) == s)
        if eqn.type == 'T'
            eqn.U( : , s + 1 : s + m) = eqn.B;
        else
            eqn.U( : , s + 1 : s + m) = eqn.C';
        end
    end
    eqn.haveUV=1;
    j=j+1;
    if restarted
        already_restarted = 1;
        restarted = 0;
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ((opts.nm.restol&&res(j-1)<opts.nm.restol)...
            ||(opts.nm.rctol&&rc(j-1)<opts.nm.rctol))
        break
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set tolerance for next ADI iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~opts.adi.inexact && (j > 2)
        % Check whether inexact ADI iteration can be turned on again
        switch opts.nm.inexact
            case 'linear'
                opts.adi.inexact = res(j - 1) < opts.nm.tau * res(j - 2);
            case 'superlinear'
                opts.adi.inexact = res(j - 1) < 1 / (j ^ 3 + 1) * res(j - 2);
            case 'quadratic'
                if res(j - 2) > 1
                    opts.adi.inexact = res(j - 1) < opts.nm.tau / sqrt(res(j - 2));
                else
                    opts.adi.inexact = res(j - 1) < opts.nm.tau * res(j - 2) * res(j - 2);
                end
        end
        if opts.adi.inexact
            warning('MESS:lrnm', ...
                'Turning inexact ADI iteration back on.');
        end
    end
    if opts.adi.inexact
        switch opts.nm.inexact
            case 'linear'
                opts.adi.outer_tol = opts.nm.tau * res(j - 1);
            case 'superlinear'
                opts.adi.outer_tol = 1 / (j ^ 3 + 1) * res(j - 1);
            case 'quadratic'
                if res(j - 1) > 1
                    opts.adi.outer_tol = opts.nm.tau / sqrt(res(j - 1));
                else
                    opts.adi.outer_tol = opts.nm.tau * res(j - 1) * res(j - 1);
                end
        end
        opts.adi.outer_tol = max(opts.adi.outer_tol, opts.adi.restol);
        opts.adi.outer_res = res(j - 1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eqn,opts,oper]=oper.mul_E_post(eqn,opts,oper);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.adi=adiout;
out.nm.niter=j - 1;
if eqn.setUV
    out.nm.K = eqn.V(:,end-m+1:end);
else
    out.nm.K = eqn.V;
end
if opts.nm.LDL_T
    out.nm.D = adiout.D;
    out.nm.S = eqn.S;
end
if opts.nm.restol, out.nm.res=res(1:out.nm.niter); end
if opts.nm.rctol, out.nm.rc=rc(1:out.nm.niter); end
