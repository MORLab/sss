function X=sol_E_so_2(eqn, opts,opE,B,opB)

% function X=sol_E_so_2(eqn, opts,opE,B,opB)
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
% This function returns X= E\B, where matrix E given by structure eqn and input matrix B could be transposed.
%
%   Inputs:
%   eqn     structure containing data for matrix E (fields 'D_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opE     character specifying the shape of E
%           opE = 'N' solves E *X = opB(B)
%           opE = 'T' sovles E'*X = opB(B)
%   B       p-x-q matrix 
%   opB     character specifying the shape of B
%           opB = 'N' solves opE(E)*X = B
%           opB = 'T' solves opE(E)*X = B'
%
%   Output:
%
%                                       |D  M|
%   X       matrix fullfilling equation |M  0| *X = opB(B)
%
%   This function does not use other so3 functions.
%
% ATTENTION: opE is not used since matrix E is symmetric

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


%% check input parameters
if (~ischar(opE) || ~ischar(opB))
    error('MESS:error_arguments', 'opE or opB is not a char');
end
opE = upper(opE); opB = upper(opB);
if(~(opE=='N' || opE=='T'))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(~(opB=='N' || opB=='T'))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if (~isnumeric(B)) || (~ismatrix(B))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(~isfield(eqn,'D_') || ~isnumeric(eqn.D_) ...
        || ~isfield(eqn,'M_') || ~isnumeric(eqn.M_))
    error('MESS:error_arguments',...
    'E consists of D and M, field eqn.D_ or eqn.M_ is not defined');
end
if ~isfield(eqn,'K_') || ~isnumeric(eqn.K_)
    error('MESS:error_arguments',...
    'Field eqn.K_ is not defined or corrupted');
end


rowK = size(eqn.K_, 1);
rowE = 2*rowK;

%% perform solve operations
switch opB
    
    % implement solve E*X = B
    case 'N'
        
        if(rowE ~= size(B,1))
            error('MESS:error_arguments','number of rows of E differ with number of rows of B')
        end
        
        X1 = eqn.M_\B(rowK+1:end,:);
        X2 = eqn.M_\(B(1:rowK,:)-eqn.D_*X1);
        X  = [X1;X2];
        
    % implement solve E*X = B'
    case 'T'
        
        if(rowE ~= size(B,2))
            error('MESS:error_arguments','number of rows of E differs with number of columns of B')
        end
        
        X1 = eqn.M_\B(:,rowK+1:end)';
        X2 = eqn.M_\(B(:,1:rowK)' - eqn.D_*X1);
        X  = [X1;X2];
        
end

end
