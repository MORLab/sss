function [Z, D] = mess_lyap(A, B, C, S, E)
% mess_lyap Solve continuous-time Lyapunov equations with 
%           sparse coefficients.
% 
%    Z = mess_lyap(A, B) solves the Lyapunov matrix equation:
% 
%        A*Z*Z' + Z*Z'*A' + B*B' = 0
% 
%    [Z, D, Y] = mess_lyap(A, B, C) solves the Sylvester equation:
% 
%        A*Z*D*Z' + Z*D*Y'*B + C = 0 (NOT YET IMPLEMENTED)
% 
%    Z = mess_lyap(A, B, [], [], E) solves the generalized Lyapunov 
%        equation:
% 
%        A*Z*Z'*E' + E*Z*Z'*A' + B*B' = 0 
% 
%    [Z, D] = mess_lyap(A, B, [], S) solves the Lyapunov matrix equation
%       in ZDZ^T formulation:
% 
%        A*Z*D*Z' + Z*D*Z'*A' + B*S*B' = 0
% 
%    [Z, D] = mess_lyap(A, B, [], S, E) solves the generalized Lyapunov 
%       equation in ZDZ^T formulation:
% 
%        A*Z*D*Z'*E' + E*Z*D*Z'*A' + B*S*B' = 0 
% 
%    For the dense case see also lyap, lyapchol, dlyap.

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

%% Usfs
oper=operatormanager('default');

%% Options
ni=nargin;
no = nargout;
if ni < 4
    S = [];
end
opts.adi.info=0;
opts.adi.restol=1e-12;
opts.adi.rctol=0;
opts.adi.maxiter=100;
opts.adi.shifts.kp=50;
opts.adi.shifts.km=25;
opts.adi.shifts.method = 'projection';
opts.adi.shifts.l0 = 6;
opts.adi.norm = 'fro';
opts.adi.computeZ = 1;


%% Equation type
eqn.type = 'N';
if no == 1
    if ~isempty(S)
        warning('MESS:ignored',...
            'Fourth argument is supposed to be empty. Data is ignored.');
    end
    eqn.A_ = A;
    eqn.G = B;
    if ni == 2
        eqn.haveE = 0;
    elseif ni == 5
        if ~isempty(C)
            warning('MESS:ignored',...
                'Third argument is supposed to be empty. Data is ignored.');
        end
        eqn.haveE = 1;
        eqn.E_ = E;
    else
        error('MESS:notimplemented', 'Feature not yet implemented!');
    end
elseif no == 2 % ZDZ^T case
    opts.adi.LDL_T = 1;
    eqn.A_ = A;
    eqn.G = B;
    eqn.S = S;
    if ni == 4
        if ~isempty(C)
            warning('MESS:ignored',...
                'Third argument is supposed to be empty. Data is ignored.');
        end
        eqn.haveE = 0;
    elseif ni == 5
        if ~isempty(C)
            warning('MESS:ignored',...
                'Third argument is supposed to be empty. Data is ignored.');
        end
        eqn.haveE = 1;
        eqn.E_ = E;
    else
        error('MESS:notimplemented', 'Feature not yet implemented!');
    end
else
    error('MESS:notimplemented', 'Feature not yet implemented!');
end
eqn.B = B;

%% Shift parameter
opts.adi.shifts.p = mess_para(eqn, opts, oper);

%% Solve Equation
[Z, out] = mess_lradi(eqn, opts, oper);

%% Prepare output
if no == 2 % ZDZ^T case 
    D = kron(diag(out.D), eqn.S);
end
