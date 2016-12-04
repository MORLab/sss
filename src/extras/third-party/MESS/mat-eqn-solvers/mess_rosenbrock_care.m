function [out] = mess_rosenbrock_care(eqn, opts, oper)
%% [out] = mess_rosenbrock_care(eqn, opts, oper)
%   LDL^T factored Rosenbrock method solving diffenrential Riccati equation
%   E*d/dt X(t)*E' =-B*B' - E*X(t)*A' - A*X(t)*E' + E*X(t)*C'*C*X(t)*E' (N)
%   E'*d/dt X(t)*E =-C'*C - E'*X(t)*A - A'*X(t)*E + E'*X(t)*B*B'*X(t)*E (T)
%   backward in time.
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
%   out     struct contains solutions for every time step
%
% Input fields in struct eqn:
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
%
%   eqn.type    possible  values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional)
%
%   eqn.haveE   possible  values: 0, 1, false, true
%               if haveE = 0: matrix E is assumed to be the identity
%               (optional)
%
%   Depending on the operator chosen by the operatormanager, additional
%   fields may be needed. For the "default", e.g., eqn.A_ and eqn.E_ hold
%   the A and E matrices. For the second order types these are given
%   implicitly by the M, D, K matrices stored in eqn.M_, eqn.D_ and eqn.K_,
%   respectively.
%
% Input fields in struct opts:
%   opts.rosenbrock.time_steps  possible  values: (N x 1) array
%                               array containing the time steps
%
%   opts.rosenbrock.stage       possible  values: 1, 2
%                               use 1-stage or 2-stage Rosenbrock method
%                               (optional)
%
%   opts.rosenbrock.info        possible  values: 0, 1, false, true
%                               turn on (1) or off (0) the status output in
%                               every Rosenbrock iteration step
%                               (optional)
%
%   opts.rosenbrock.gamma       possible  values: scalar
%                               Rosenbrock method coefficient
%                               (optional)
%                               
%
%   opts.rosenbrock.save_solution
%                               possible  values: 0, 1, false, true
%                               save only K (0) or also the solution
%                               factors L and D (1)
%                               (optional)
%                               
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
%                               every iteration ADI step
%                               (optional)
%
%   opts.adi.norm               possible  values: 2, 'fro'
%                               use 2-norm (2) or Frobenius norm ('fro') to
%                               compute residual and relative change norms
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
% Output fields in struct out:
%
%   out.Ks  cell array with matrix K for every time step
%
%   out.Ls  cell array with solution factor L for every time step
%           (only if opts.rosenbrock.save_solution = 1)
%
%   out.Ds  cell array with solution factor D for every time step
%           (only if opts.rosenbrock.save_solution = 1)                             
%
% If optional input arguments are missing they may be set to default values
% and a 'MESS:control_data' warning is printed. to turn warnings off use
% warning('OFF', 'MESS:control_data')
%
%
%   Note: Uses mess_lradi to solve inner Lyapunov Equations.
%   
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
% Check for Rosenbrock Control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(opts,'rosenbrock') || ~isstruct(opts.rosenbrock)
    error('MESS:control_data','Rosenbrock control structure opts.rosenbrock missing.');
end % Single fields are checked below or inside mess_lradi
if ~isfield(opts.rosenbrock,'time_steps')
    error('MESS:control_data','opts.rosenbrock.time_steps is missing.');
end
if ~isfield(opts.rosenbrock,'stage'),  opts.rosenbrock.stage=1; end
if (opts.rosenbrock.stage ~= 1) && (opts.rosenbrock.stage ~= 2)
    error('MESS:control_data',['opts.rosenbrock.stage has to be 1 or 2.',...
        ' Other stages are not implemented']);
end
if ~isfield(opts.rosenbrock,'info'),  opts.rosenbrock.info=0; end
if ~isfield(opts.rosenbrock,'column_compression_tol'),  ...
        opts.rosenbrock.column_compression_tol = eps * oper.size(eqn, opts); end
if ~isfield(opts.rosenbrock,'gamma'),  ...
        opts.rosenbrock.gamma=1 + 1 / sqrt(2); end
if ~isfield(opts.rosenbrock,'save_solution')  
    opts.rosenbrock.save_solution = 0; 
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for ADI control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts,'adi') || ~isstruct(opts.adi)
    error('MESS:control_data','ADI control structure opts.adi missing.');
end
if ~isfield(opts.adi,'computeZ')||~isnumeric(opts.adi.computeZ) || ...
        (~opts.adi.computeZ)
    warning('MESS:control_data', ...
        'Missing or Corrupted computeZ field. Switching to default.');
    opts.adi.computeZ = 1;
end
if ~isfield(opts.adi,'accumulateK')||~isnumeric(opts.adi.accumulateK) || ...
        (~opts.adi.accumulateK) || ~opts.adi.accumulateK
    warning('MESS:control_data', ...
        'Missing or Corrupted accumulateK field. Switching to default.');
    opts.adi.accumulateK = 1;
end
if ~isfield(opts.adi.shifts,'period')
    opts.adi.shifts.period=1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(eqn, 'C') || ~isnumeric(eqn.C)
    error('MESS:control_data', 'eqn.C is not defined or corrupted');
end
if ~isfield(eqn, 'B') || ~isnumeric(eqn.B)
    error('MESS:control_data', 'eqn.B is not defined or corrupted');
end
if ~isfield(eqn, 'L0') || ~isnumeric(eqn.L0)
    eqn.L0 = [];
end
if ~isfield(eqn, 'type')
    eqn.type = 'N';
    warning('MESS:MESS:control_data',['Unable to determine type of equation.'...
        'Falling back to type ''N''']);
elseif (eqn.type ~= 'N') && (eqn.type ~= 'T')
    error('MESS:equation_type', 'Equation type must be either ''T'' or ''N''');
end
[eqn, erg] = oper.init(eqn, opts, 'A','E');
if ~erg
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end
if eqn.type == 'T'
    q = size(eqn.C, 1);
else
    q = size(eqn.B, 2);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All checks done. Here comes the real work!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flip time_steps to run backwards in time
times = opts.rosenbrock.time_steps(end:-1:1);

L = eqn.L0;
D = eye(size(L, 2));
if ~isempty(L)
    if eqn.type == 'T'
        K = oper.mul_E(eqn, opts, eqn.type, (L * (L' * eqn.B)), 'N');
    else
        K = oper.mul_E(eqn, opts, eqn.type, (L * (L' * eqn.C)), 'T');
    end
else
    K = 0;
    BLD = [];
end
Iq = eye(q);
if ~isempty(L)
    if opts.rosenbrock.save_solution
        out.Ls = {L}; % L of step 0
        out.Ds = {D}; % D of step 0
    end
    out.Ks = {K}; % K of step 0
else
    if opts.rosenbrock.save_solution
        out.Ls = {0}; % L of step 0
        out.Ds = {0}; % D of step 0
    end
    out.Ks = {0};
end
if eqn.type == 'T'
    eqn.G = eqn.C';
else
    eqn.G = eqn.B;
end
G = eqn.G;
eqn.haveUV = 0;
opts.adi.LDL_T = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2 : length(times)
    t = times(j);
    opts.rosenbrock.tau = times(j - 1) - times(j);
    if opts.rosenbrock.stage == 1
        if eqn.type == 'T'
            eqn.U = eqn.B;
        else
            eqn.U = eqn.C';
        end
    else
        if eqn.type == 'T'
            eqn.U = opts.rosenbrock.tau * opts.rosenbrock.gamma * eqn.B;
        else
            eqn.U = opts.rosenbrock.tau * opts.rosenbrock.gamma * eqn.C';
        end
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update E^T * L
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(L)
        EL = oper.mul_E(eqn, opts, eqn.type, L, 'N');
        m = size(EL, 2);
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update DLB block for S
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if eqn.type == 'T'
            BLD_tmp = (eqn.B' * L) * D;
        else
            BLD_tmp = (eqn.C * L) * D;
        end
        if opts.rosenbrock.stage == 1
            BLD = BLD_tmp' * BLD_tmp + 1/opts.rosenbrock.tau * D;
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update G and S
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            eqn.G = [G, EL];
        else % p = 2
            BLD = [zeros(m ,m), D;
                D, -BLD_tmp' * BLD_tmp];
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update G and S
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            eqn.G =  [G, oper.mul_A(eqn, opts, eqn.type, L, 'N'), EL];
        end
    end
    eqn.S = blkdiag(Iq, BLD);
    [eqn.G, eqn.S] = mess_LDL_column_compression(eqn.G, 'N',eqn.S, 'N',...
            opts.rosenbrock.column_compression_tol);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute new ADI shifts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~mod(j - 2,opts.adi.shifts.period)
        opts.adi.shifts.p=mess_para(eqn,opts,oper);
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform the actual step computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.rosenbrock.stage == 1
        [L, adiout]=mess_lradi(eqn,opts, oper);
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform column compression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [L, D] = mess_LDL_column_compression(L, 'N', kron(diag(adiout.D), eqn.S), ...
            'N', opts.rosenbrock.column_compression_tol);
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract K
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K = adiout.Knew;
    else % stage = 2
        [T1, adiout1]=mess_lradi(eqn,opts, oper);
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform column compression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [T1, D1] = mess_LDL_column_compression(T1, 'N', kron(diag(adiout1.D), eqn.S),...
            'N', opts.rosenbrock.column_compression_tol);
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update RHS for second equation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ET1 = oper.mul_E(eqn, opts, eqn.type, T1, 'N');
        eqn.G = ET1;
        if eqn.type == 'T'
            BT1D_tmp = (eqn.B' * T1) * D1;
        else
            BT1D_tmp = (eqn.C * T1) * D1;
        end
        BT1D = opts.rosenbrock.tau * opts.rosenbrock.tau * ...
            BT1D_tmp' * BT1D_tmp + ...
            (2 - 1/opts.rosenbrock.gamma) * D1;
        eqn.S = BT1D;
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute new ADI shifts
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(j - 2,opts.adi.shifts.period)
            opts.adi.shifts.p=mess_para(eqn,opts,oper);
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % solve second equation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [T2, adiout2]=mess_lradi(eqn,opts, oper);
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct new X = L * D * L^T
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L = [L, T1, T2];
        D = blkdiag(D, opts.rosenbrock.tau * (2 - 1 / (2 * opts.rosenbrock.gamma)) * D1, ...
            - opts.rosenbrock.tau / 2. * kron(diag(adiout2.D), eqn.S));
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform column compression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [L, D] = mess_LDL_column_compression(L, 'N', D, 'N', ...
            opts.rosenbrock.column_compression_tol);
        if eqn.type == 'T'
            eqn.G = eqn.C';
        else
            eqn.G = eqn.B;
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct new K
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if opts.adi.accumulateK
            K = K + opts.rosenbrock.tau * (2 - 1 / (2 * opts.rosenbrock.gamma)) * adiout1.Knew ...
                - opts.rosenbrock.tau / 2 * adiout2.Knew;
        else
            if eqn.type == 'T'
                K = oper.mul_E(eqn, opts,eqn.type,L,'N') * (D * (L' * eqn.B));
            else
                K = oper.mul_E(eqn, opts,eqn.type,L,'N') * (D * (L' * eqn.C'));
            end
        end
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print status information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.rosenbrock.info
        fprintf('rosenbrock step %4d s\n', t);
    end
    if opts.rosenbrock.save_solution
        out.Ls = [{L}, out.Ls];
        out.Ds = [{D}, out.Ds];
    end
    out.Ks = [{K}, out.Ks];
    eqn.V = K;
    eqn.haveUV = 1;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
