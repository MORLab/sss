function [out] = mess_bdf_care(eqn, opts, oper)
%% function [out] = mess_bdf_care(eqn, opts, oper)
%   LDL^T factored BDF method solving diffenrential Riccati equation
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
%   eqn.K0      dense (n x m1) matrix, initial K
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
%   opts.bdf.time_steps         possible  values: (N x 1) array
%                               array containing the time steps
%
%   opts.bdf.step               possible  values: 1, 2, 3, 4, 5, 6
%                               use p-step  Rosenbrock method with 1<=p<=6
%                               (optional)
%
%   opts.bdf.info               possible  values: 0, 1, false, true
%                               turn on (1) or off (0) the status output in
%                               every Rosenbrock iteration step
%                               (optional)
%
%   opts.bdf.save_solution      possible  values: 0, 1, false, true
%                               save only K (0) or also the solution
%                               factors L and D (1)
%                               (optional)
%                               
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
% Output fields in struct out:
%
%   out.Ks  cell array with matrix K for every time step
%
%   out.Ls  cell array with solution factor L for every time step
%           (only if opts.bdf.save_solution = 1)
%
%   out.Ds  cell array with solution factor D for every time step
%           (only if opts.bdf.save_solution = 1)
%                               
%
% If optional input arguments are missing they may be set to default values
% and a 'MESS:control_data' warning is printed. to turn warnings off use
% warning('OFF', 'MESS:control_data')
%
%   Note: uses mess_lrnm to solve inner Riccati Equations
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
% Check for BDF Control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(opts,'bdf') || ~isstruct(opts.bdf)
    error('MESS:control_data','BDF control structure opts.bdf missing.');
end % Single fields are checked below or inside mess_lradi
if ~isfield(opts.bdf,'time_steps')
    error('MESS:control_data','opts.bdf.time_steps is missing.');
end

if ~isfield(opts.bdf,'step'),  opts.bdf.step=1; end
if rem(opts.bdf.step, 1) || (opts.bdf.step < 1) || (opts.bdf.step > 6)
    error('MESS:control_data','opts.bdf.step has an invalid value.');
end
if ~isfield(opts.bdf,'info'),  opts.bdf.info=0; end
if ~isfield(opts.bdf,'column_compression_tol'),  ...
    opts.bdf.column_compression_tol=eps * oper.size(eqn, opts); end
if ~isfield(opts.bdf,'save_solution'),  opts.bdf.save_solution=0; end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for Newton control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts,'nm') || ~isstruct(opts.nm)
    error('MESS:control_data','Newton control structure opts.nm missing.');
end   
opts.nm.LDL_T = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for ADI control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts,'adi') || ~isstruct(opts.adi)
    error('MESS:control_data','ADI control structure opts.adi missing.');
end 
if ~isfield(opts.adi,'computeZ')||~isnumeric(opts.adi.computeZ) || ...
    (~opts.adi.computeZ)
    warning('MESS:computeZ', ...
        'Missing or Corrupted computeZ field. Switching to default.');
    opts.adi.computeZ = 1;
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
if ~isfield(eqn, 'C') || ~isnumeric(eqn.C)
    error('MESS:control_data', 'eqn.C is not defined or corrupted');
end
if ~isfield(eqn, 'B') || ~isnumeric(eqn.B)
    error('MESS:control_data', 'eqn.B is not defined or corrupted');
end
if ~isfield(eqn, 'L0') || ~isnumeric(eqn.L0)
    eqn.L0 = [];
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
% Initialize BDF parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = zeros(6, 6);
alpha(1, 1) = -1;
alpha(1 : 2, 2) = [-4. / 3.; 1. / 3.];
alpha(1 : 3, 3) = [-18. / 11.; 9. / 11.; -2. / 11.];
alpha(1 : 4, 4) = [-48. / 25.; 36. / 25.; -16. / 25.; 3. / 25.];
alpha(1 : 5, 5) = [-300. / 137.; 300. / 137.; -200. / 137.; 75. / 137.; ...
    -12. / 137.];
alpha(1 : 6, 6) = [-360. / 147.; 450. / 147.; -400. / 147.; 225. / 147.; ...
    -72. / 147.; 10. / 147.];
beta = [1; 2 / 3; 6 / 11; 12 / 25; 60 / 137; 60 / 147];
step = opts.bdf.step;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
times = opts.bdf.time_steps(end:-1:1);
L = eqn.L0;
D = eye(size(L, 2));
if ~isempty(L)
    if eqn.type == 'T'
        K0 = oper.mul_E(eqn, opts, eqn.type, (L * (L' * eqn.B)), 'N');
    else
        K0 = oper.mul_E(eqn, opts, eqn.type, (L * (L' * eqn.C')), 'N');
    end
end

L_lenghts = zeros(opts.bdf.step, 1);
Ds = {};
alphaDs = [];
Iq = eye(q);
if ~isempty(L)
    if opts.bdf.save_solution
        out.Ls = {L}; % L of step 0
        out.Ds = {D}; % D of step 0
    end
    out.Ks = {K0}; % K of step 0
else
    if opts.bdf.save_solution
        out.Ls = {0}; % L of step 0
        out.Ds = {0}; % D of step 0
    end
    out.Ks = {0}; % K of step 0
end
if eqn.type == 'T'
    C = eqn.C;
else
    B = eqn.B;
end
ETL = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2 : length(times)
    t = times(j);
    opts.bdf.tau = times(j - 1) - times(j);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set order and beta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opts.bdf.step = min(step, j - 1);
    opts.bdf.beta = beta(opts.bdf.step);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update E' * L
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(L)
        if isempty(ETL)
            ETL = oper.mul_E(eqn, opts, eqn.type, L, 'N');
            L_lenghts(1) = size(L, 2);
        else
            ETL = [oper.mul_E(eqn, opts, eqn.type, L, 'N'), ...
                ETL( : , 1 : end - L_lenghts(step))];
            L_lenghts(2 : step) = L_lenghts(1 : end - 1);
            L_lenghts(1) = size(L, 2);
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update D blocks for S
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alphaDs = -alpha(1, opts.bdf.step) * D;
        for i = 2 : opts.bdf.step
            if size(Ds, 2) >= i - 1
                alphaDs = blkdiag(alphaDs, ...
                    -alpha(i, opts.bdf.step) * Ds{i - 1});
            end
        end
        if size(Ds, 2) == opts.bdf.step
            Ds = {D, Ds{1 : end - 1}};
        else
            Ds = [{D}, Ds];
        end
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update C/B and S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if eqn.type == 'T'
        eqn.C = [C; ETL'];
    else
        eqn.B = [B, ETL];
    end
    eqn.S = blkdiag(opts.bdf.tau * opts.bdf.beta * Iq, alphaDs);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform column compression
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if eqn.type == 'T'
        [eqn.C, eqn.S] = mess_LDL_column_compression(eqn.C, 'T', eqn.S,...
            'N', opts.bdf.column_compression_tol);
    else
        [eqn.B, eqn.S] = mess_LDL_column_compression(eqn.B, 'N', eqn.S,...
            'N', opts.bdf.column_compression_tol);
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform the actual step computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [L, nmout]=mess_lrnm(eqn,opts, oper);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform column compression
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [L, D] = mess_LDL_column_compression(L, 'N', kron(diag(nmout.nm.D), ...
        nmout.nm.S), 'N', opts.bdf.column_compression_tol);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print status information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.bdf.info
        fprintf('BDF step %4d s\n', t);
        fprintf('\t Newton steps: %2d \n', nmout.nm.niter);
    end
    if opts.bdf.save_solution
        out.Ls = [{L}, out.Ls];
        out.Ds = [{D}, out.Ds];
    end
    out.Ks = [{nmout.nm.K}, out.Ks];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
