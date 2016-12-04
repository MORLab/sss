function n = size_so_2(eqn, opts)

% function n = size_so_2(eqn, opts)
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
%% check data

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
if ~isfield(eqn,'K_') || ~isnumeric(eqn.K_)
    error('MESS:error_arguments',...
    'Field eqn.K_ is not defined or corrupted');
end

%% compute size
n = 2*size(eqn.K_,1);

end
