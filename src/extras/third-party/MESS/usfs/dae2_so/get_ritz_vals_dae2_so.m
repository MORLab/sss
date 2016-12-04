function [rw, Hp, Hm, Vp, Vm] = get_ritz_vals_dae2_so(eqn, opts, oper, U, W, p_old)
%  This function is an exact copy of the Penzl heuristic part in mess_para.
%  the only difference is that B or C and K are filled up by trailing zero
%  blocks to allow for the computation of the Ritz values with respect to
%  the full size block matrices instead of the restriction to the
%  (1,1)-block for which the ADI is formulated.
%
%  MMESS (Jens Saak, October 2013)

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

% Input data not completely checked!
for mat='MDKG'
    if(~isfield(eqn, sprintf('%c_',mat)) || ~eval(sprintf('isnumeric(eqn.%c_)',mat)))
        error('MESS:error_arguments', 'field eqn.%c_ is not defined',mat);
    end
end

[eqn, erg] = oper.init(eqn, opts, 'A','E');
if ~erg
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end
% returns order of A or states of A, A is supposed to be square
nv = size(eqn.M_,1);
np = size(eqn.G_,1);

%%
% here we add the trailing zero blocks. Note that we are not passing the
% eqn structure back as an output, so this change is not visible in
% anything above this routine and will only be passed on to the function
% handles used in here.
if isfield(eqn, 'U') && ~isempty(eqn.U)
    eqn.U = [eqn.U; sparse(np, size(eqn.U, 2))];
end
if isfield(eqn,'V') && ~isempty(eqn.V)
    eqn.V=[eqn.V; sparse(np, size(eqn.V,2))]; 
end

if (~isfield(opts.adi.shifts, 'b0') || isempty(opts.adi.shifts.b0))
    opts.adi.shifts.b0 = ones(n,1);
else 
    if length(opts.adi.shifts.b0) ~= 2*nv+np
        warning('MESS:b0',...
        'b0 has the wrong length. Switching to default.');
        opts.adi.shifts.b0 = ones(2*nv+np,1);
    end
end
if isfield(opts.adi.shifts, 'method') && ...
        strcmp(opts.adi.shifts.method, 'projection')
    U = [U; zeros(2 * nv + np - size(U, 1), size(U, 2))];
    if isempty(W)
        % first shifts are computed with U = eqn.G and W = A * eqn.G
        W = oper.mul_A(eqn, opts, eqn.type, U, 'N');
    else
        W = [W; zeros(2 * nv + np - size(W, 1), size(W, 2))];
    end
    rw = mess_projection_shifts(eqn, opts, oper, U, ...
        W, p_old);
else
    [rw, Hp, Hm, Vp, Vm] = mess_get_ritz_vals(eqn, opts, oper);
end
end
