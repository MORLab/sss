function n = size_so_1(eqn, opts)

% function n = size_so_1(eqn, opts)
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
% This function returns the number of rows of the matrices A and E.
%
%    Inputs:
%
%    eqn       structure containing field 'K_'
%    opts      structure containing parameters for the algorithm
%
%    Output:
%
%    n         double number of rows of matrix V_ in structure eqn
%
% This function does not use other so1 functions.

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
if(~isfield(eqn,'K_') || ~isnumeric(eqn.K_))
    error('MESS:error_arguments',...
        'A consists of K and D, field eqn.K_ is not defined or corrupted');
end
n = 2*size(eqn.K_,1);

end
