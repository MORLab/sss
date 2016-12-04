function [ W, res0 ] = init_res_so_1( eqn, opts, RHS)

% function [ W, res0 ] = init_res_so_1( eqn, opts, RHS)
%
% The second order system
%
%    M x'' + D x' + K x = B u
%                       y = C x
%
% is transformed to the first order system 
%
%    E x' = A x + B u
%
% where
% 
%       |-K  0 |
%    E= | 0  M | ,
%
%       | 0 -K |
%    A= |-K -D |,
%
%       | 0 |
%    B= | B |,
%
%       | x |
%    x= | x'|.
%
% Matrices M, D, K are assumed to be symmetric and quadratic.
% Matrix K has full rank.
%
%
% This function returns a matrix W and its associated residuum res0.
%
%   Inputs:
%
%   eqn        structure containing data for G or B or C
%   opts       structure containing parameters for the algorithm
%   RHS        right hand side matrix
% 
%   Outputs:
% 
%   W          matrix given by ADI to compute residuum
%   res0       initial residuum
%
% This function does not use other so1 functions.
%% check input data

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
if (~isnumeric(RHS)) || (~ismatrix(RHS))
    error('MESS:error_arguments','RHS has to ba a matrix');
end
%% compute low rank residual
W = RHS;
    
%% compute res0
if opts.adi.LDL_T
    res0 = max(abs(eig(W' * W * eqn.S)));
else
    res0 = norm(W' * W, 2);
end

end

