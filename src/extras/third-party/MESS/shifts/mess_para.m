function [p, err_code, rw, Hp, Hm, Vp, Vm, eqn, opts, oper] = mess_para(eqn, opts, oper)
%
%  Estimation of suboptimal ADI shift parameters for the matrix (operator) F=A
%
%  Calling sequence:
%
%    [p, err_code, rw, Hp, Hm, Vp, Vm, eqn, opts, oper] = mess_para(eqn, opts, oper)
%
%  Input:
%
%    eqn       structure contains data A, E, B, C, K
%
%    opts      struct contains parameters for the algorithm
%
%    oper      contains function handles with operations for A and E
%
%  Output:
%
%    p         an opts.adi.shifts.l0- or opts.adi.shifts.l0+1-vector of
%              suboptimal ADI parameters;
%    err_code  Error code; = 1, if Ritz values with positive real parts
%              have been encountered; otherwise, err_code = 0;
%    rw        vector containing the Ritz values;
%    Hp        Hessenberg matrix in Arnoldi process w.r.t. F;
%    Hm        Hessenberg matrix in Arnoldi process w.r.t. inv(F);
%    Vp        Orthogonal matrix in Arnoldi process w.r.t. F;
%    Vm        Orthogonal matrix in Arnoldi process w.r.t. inv(F);
%
% Input fields in struct eqn:
%
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
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
%   opts.adi.shifts.l0          possible  values: integer > 0
%                               number of shifts that should be computed
%                               2*l0 < kp + km is required
%                               (optional)
%
%   opts.adi.shifts.kp          possible  values: integer > 0
%                               number of Arnoldi steps w.r.t. F for
%                               heuristic shift computation
%                               kp < n is required
%                               (optional)
%
%   opts.adi.shifts.km          possible  values: integer > 0
%                               number of Arnoldi steps w.r.t. inv(F) for
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
%  Remarks:
%
%    Typical values are opts.adi.shifts.l0 = 10..40,
%    opts.adi.shifts.kp = 20..80, opts.adi.shifts.km = 10..40.
%    The harder the problem is the large values are necessary.
%    Larger values mostly result in a faster convergence, but also in a
%    larger memory requirement.
%    However, for "well-conditioned" problems small values of
%    opts.adi.shifts.l0 can lead to the optimal performance.
%
%  References:
%
%  [1] T. Penzl.
%      LYAPACK (Users' Guide - Version 1.0).
%      1999.
%
%   uses operatorfunctions size directly and indirectly
%   size, sol_A, mul_A, sol_E, mul_E in mess_arn

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

%  MMESS (Jens Saak, October 2013)

% Input data not completely checked!

%% check data
if ~isfield(opts,'adi') || ~isstruct(opts.adi)
    error('MESS:control_data','ADI control structure opts.ADI missing.');
end
if ~isfield(opts.adi,'shifts') || ~isstruct(opts.adi.shifts)
    warning('MESS:control_data',['shift parameter control structure missing.', ...
        'Switching to default l0 = 25, kp = 50, km = 25.']);
    opts.adi.shifts.l0 = 25;
    opts.adi.shifts.kp = 50;
    opts.adi.shifts.km = 25;
else
    if ~isfield(opts.adi.shifts,'l0')||~isnumeric(opts.adi.shifts.l0)
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.adi.shifts.l0 field.', ...
            'Switching to default: 25']);
        opts.adi.shifts.l0 = 25;
    end
    if ~isfield(opts.adi.shifts,'kp')||~isnumeric(opts.adi.shifts.kp)
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.adi.shifts.kp field.', ...
            'Switching to default: 50']);
        opts.adi.shifts.kp = 50;
    end
    if ~isfield(opts.adi.shifts,'km')||~isnumeric(opts.adi.shifts.km)
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.adi.shifts.km field.', ...
            'Switching to default: 25']);
        opts.adi.shifts.km = 25;
    end
end
if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
if ~isfield(eqn, 'type'), eqn.type = 'N'; end
[eqn, erg] = oper.init(eqn, opts, 'A','E');
if ~erg
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end
err_code = 0;
if ~isfield(opts,'rosenbrock'), opts.rosenbrock=[]; end
if isstruct(opts.rosenbrock)&&isfield(opts.rosenbrock,'tau')
    rosenbrock = 1;
    if opts.rosenbrock.stage == 1
        pc = -1 / (2 * opts.rosenbrock.tau);
        taugamma = 1;
    else % p = 2
        taugamma = (opts.rosenbrock.tau * opts.rosenbrock.gamma);
        pc = ( - 0.5) / taugamma;
    end
else
    rosenbrock = 0;
end
if ~isfield(opts,'bdf'), opts.bdf=[]; end
if isstruct(opts.bdf) && isfield(opts.bdf, 'tau') && isfield(opts.bdf, 'beta')
    bdf = 1;
    pc = -1 / (2 * opts.bdf.tau * opts.bdf.beta);
else
    bdf = 0;
end

%% initialize usfs
[eqn,opts, oper] = oper.mul_A_pre(eqn, opts, oper);
[eqn,opts, oper] = oper.mul_E_pre(eqn, opts, oper);
[eqn,opts, oper] = oper.sol_A_pre(eqn, opts, oper);
[eqn,opts, oper] = oper.sol_E_pre(eqn, opts, oper);
%%
if ~isfield(opts.adi.shifts,'method')
    opts.adi.shifts.method='heur';
end

switch opts.adi.shifts.method
    case {'heur','heuristic','penzl','Penzl'}
        %%
        if isfield(oper,'get_ritz_vals')
            if nargout < 4
                rw = oper.get_ritz_vals(eqn,opts,oper);
            else
                [rw, Hp, Hm, Vp, Vm, eqn, opts, oper] = oper.get_ritz_vals(eqn,opts,oper);
            end
        else
            if nargout < 4
                rw = mess_get_ritz_vals(eqn,opts,oper);
            else
                [rw,  Hp, Hm, Vp, Vm, eqn, opts, oper] = mess_get_ritz_vals(eqn,opts,oper);
            end
        end
        
        p = mess_mnmx(rw,opts.adi.shifts.l0);
        
    case {'wachspress','Wachspress'}
        %%
        if isfield(oper,'get_ritz_vals')
            if nargout < 4
                rw = oper.get_ritz_vals(eqn,opts,oper);
            else
                [rw, Hp, Hm, Vp, Vm, eqn, opts, oper] = oper.get_ritz_vals(eqn,opts,oper);
            end
        else
            if nargout < 4
                rw = mess_get_ritz_vals(eqn,opts,oper);
            else
                [rw,  Hp, Hm, Vp, Vm, eqn, opts, oper] = mess_get_ritz_vals(eqn,opts,oper);
            end
        end
        
        a=min(abs(real(rw)));
        b=max(abs(real(rw)));
        alpha=atan(max(imag(rw)./real(rw)));
        
        if ~isfield(opts.adi.shifts,'wachspress')
            opts.adi.shifts.wachspress='T';
        end
        switch opts.adi.shifts.wachspress
            case 'N'
                p = mess_wachspress_n(a,b,alpha,opts.adi.shifts.l0);
            case 'T'
                p = mess_wachspress(a,b,alpha,opts.adi.restol);
            otherwise
                error('MESS:shift_method','wachspress selector needs to be either ''T'' or ''N''');
        end
    case 'projection'
        if nargout > 3
            error('MESS:shift_method', ...
                'For shift method ''projection'' matrices Hp, Hm, Vp and Vm are not available.')
        end
        if isfield(eqn, 'G')
            U = eqn.G;
        elseif eqn.type == 'N'
            U = eqn.B;
        else
            U = eqn.C';
        end
        if issparse(U)
            U = full(U);
        end
        p = [];
        i = 1;
        while isempty(p)
            if bdf
                if isfield(oper,'get_ritz_vals')
                    p = oper.get_ritz_vals(eqn, opts, oper, U, ...
                        (opts.bdf.tau * opts.bdf.beta) * ...
                        oper.mul_ApE(eqn, opts, eqn.type, pc, eqn.type, U, 'N'), []);
                else
                    p = mess_projection_shifts(eqn, opts, oper, U, ...
                        (opts.bdf.tau * opts.bdf.beta) * ...
                        oper.mul_ApE(eqn, opts, eqn.type, pc, eqn.type, U, 'N'), []);
                end
            elseif rosenbrock
                if eqn.haveUV
                    if eqn.type == 'N'
                        if isfield(oper,'get_ritz_vals')
                            p = oper.get_ritz_vals(eqn, opts, oper, U, ...
                                taugamma * oper.mul_ApE(eqn, opts, eqn.type, pc, eqn.type, U, 'N') ...
                                - eqn.U * (eqn.V' * U) , []);
                        else
                            p = mess_projection_shifts(eqn, opts, oper, U, ...
                                taugamma * oper.mul_ApE(eqn, opts, eqn.type, pc, eqn.type, U, 'N') ...
                                - eqn.U * (eqn.V' * U) , []);
                        end
                    else
                        if isfield(oper,'get_ritz_vals')
                            p = oper.get_ritz_vals(eqn, opts, oper, U, ...
                                taugamma * oper.mul_ApE(eqn, opts, eqn.type, pc, eqn.type, U, 'N') ...
                                - eqn.V * (eqn.U' * U) , []);
                        else
                            p = mess_projection_shifts(eqn, opts, oper, U, ...
                                taugamma * oper.mul_ApE(eqn, opts, eqn.type, pc, eqn.type, U, 'N') ...
                                - eqn.V * (eqn.U' * U) , []);
                        end
                    end
                else
                    if isfield(oper,'get_ritz_vals')
                        p = oper.get_ritz_vals(eqn, opts, oper, U, ...
                            taugamma * oper.mul_ApE(eqn, opts, eqn.type, pc, eqn.type, U, 'N') ...
                            , []);
                    else
                        p = mess_projection_shifts(eqn, opts, oper, U, ...
                            taugamma * oper.mul_ApE(eqn, opts, eqn.type, pc, eqn.type, U, 'N') ...
                            , []);
                    end
                end
            else
                if isfield(oper,'get_ritz_vals')
                    p = oper.get_ritz_vals(eqn, opts, oper, U, ...
                        oper.mul_A(eqn, opts, eqn.type, U, 'N'), []);
                else
                    p = mess_projection_shifts(eqn, opts, oper, U, ...
                        oper.mul_A(eqn, opts, eqn.type, U, 'N'), []);
                end
            end
            if isempty(p)
                if  (i < 5)
                    warning('MESS:mess_para',['Could not compute initial projection shifts. ',...
                        'Retry with random RHS.'])
                    U = rand(size(U));
                else
                    error('MESS:mess_para','Could not compute initial projection shifts.')
                end
            end
            i = i + 1;
        end
    otherwise
        error('MESS:shift_method','unknown shift computation method requested.')
end

p = cplxpair( p, 1000*eps(p(1)) ); % ensure that complex pairs are
                                % actually paired. The tolerance is
                                % increased by a factor of 10
                                % compared to the default to ensure
                                % this also works in Octave where
                                % eig seems to be less accurate. 

%% finalize usfs
[eqn,opts, oper] = oper.mul_A_post(eqn, opts, oper);
[eqn,opts, oper] = oper.mul_E_post(eqn, opts, oper);
[eqn,opts, oper] = oper.sol_A_post(eqn, opts, oper);
[eqn,opts, oper] = oper.sol_E_post(eqn, opts, oper);
end
