function [rw,  Hp, Hm, Vp, Vm, eqn, opts, oper] = get_ritz_vals_so_2(eqn, opts, oper, U, W, p_old)

% function [rw,  Hp, Hm, Vp, Vm, eqn, opts, oper] = get_ritz_vals_so_2(eqn,opts,oper)
%
% The second order system
%
%   M x"(t) + D x'(t) + K x(t) = B u(t)
%                         y(t) = C x(t)
%       
% is transformed to the first order system
%
%   E z'(t) = A z(t) + G u(t)
%      y(t) = L z(t)
%  
% where
%
%      | D  M|
%   E= | M  0|
%   
%      |-K  0|
%   A= | 0  M|
%   
%      | B |
%   G= | 0 |
%   
%   L= [C  0]
%   
%         | x(t)  |
%   z(t)= | x'(t) | .
%   
% Matrices M, D, K are assumed to be quadratic, symmetric and positive definit.
%
% 
% This function returns suitable Ritz values, Hessenberg matrices and matrices consisting of basis vectors corresponding to A and A^{-1}. 
%
%   Input:
%
%   eqn      data structure
%   opts     structure containing parameters for the algorithm
%   oper     
%
%   Output:
%
%   rw       vector of Ritz values
%   Hp       Hessenberg matrix corresponding to A
%   Hm       Hessenberg matrix corresponding to A^{-1}
%   Vp       matrix consisting of basis vectors corresponding to A
%   Vm       matrix consisting of basis vectors corresponding to A^{-1}
%
% This function does not use other so3 functions.
% This function uses another help function mess_get_ritz_vals(eqn,opts,oper), which returns all Ritz values and Hessenberg matrices
% and matrices consisting of basis vectors corresponding to A and A^{-1}. 

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

if isfield(opts.adi.shifts, 'method') && ...
        strcmp(opts.adi.shifts.method, 'projection')
    if isempty(W)
        % first shifts are computed with U = eqn.G and W = A * eqn.G
        W = oper.mul_A(eqn, opts, eqn.type, U, 'N');
    end
    rw = mess_projection_shifts(eqn, opts, oper, U, ...
        W, p_old);
else
    [rw, Hp, Hm, Vp, Vm] = mess_get_ritz_vals(eqn, opts, oper);
end
rw = rw(abs(rw)<1e6);
rw = rw(abs(rw)>1e-6);
end
