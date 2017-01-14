function X=sol_ApE_so_2(eqn, opts,opA,p,opE,B,opB)

% function X=sol_ApE_so_2(eqn, opts,opA,p,opE,C,opC)
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
% This function returns X =(A+p*E)\C, where matrices A and E given by structure eqn and input matrix C could be transposed.
% Matrices A and E are assumed to be quadratic.
%
%   Inputs:
%
%   eqn     structure containing data for matrices A (fields 'K_' and 'M_') and E (fields 'D_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' solves (A + p* opE(E))*X = opC(C) 
%           opA = 'T' solves (A' + p* opE(E))*X = opC(C) 
%   p       scalar value
%   opE     character specifying the shape of E
%           opE = 'N' solves (opA(A) + p* E)*X = opC(C) 
%           opE = 'T' solves (opA(A) + p* E')*X = opC(C) 
%   B       n-x-p matrix
%   opB     character specifies the shape of B
%           opC = 'N' solves (opA(A) + p* opE(E))*X = B
%           opC = 'T' solves (opA(A) + p* opE(E))*X = B'
%
%   Output:
%
%                                       (|-K  0|     |D  M|)
%   X       matrix fullfilling equation (| 0  M| + p*|M  0|)*X = opB(B)
%
%   This function does not use other so3 functions.
%
% ATTENTION: opA and opE are not used since matrices A and E are symmetric

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
if (~ischar(opA) || ~ischar(opE) || ~ischar(opB))
    error('MESS:error_arguments', 'opA, opE or opC is not a char');
end

opA = upper(opA); opE = upper(opE); opB = upper(opB);

if(~(opA=='N' || opA=='T'))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(~(opE=='N' || opE=='T'))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(~(opB=='N' || opB=='T'))
    error('MESS:error_arguments','opC is not ''N'' or ''T''');
end

if(~isnumeric(p))
    error('MESS:error_arguments','p is not numeric');
end

if (~isnumeric(B)) || (~ismatrix(B))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
if(eqn.haveE ==1)
    if(~isfield(eqn,'M_') || ~isnumeric(eqn.M_) || ~isfield(eqn,'D_') || ...
            ~isnumeric(eqn.D_) || ~isfield(eqn,'K_') || ~isnumeric(eqn.K_))
        error('MESS:error_arguments',...
            'field eqn.M_, eqn.D_ or eqn.K_ is not defined or corrupted');
    end
else
    error('MESS:error_arguments','eqn.haveE has to be 1 because of the structure of E')
end

rowK = size(eqn.K_, 1);
rowA = 2*rowK;

%% perform solve operations
switch opB
    
    % implement solve (A+p*E)*X = B
    case 'N'
        
        if (rowA ~= size(B,1))
            error('MESS:error_arguments','number of rows of A differs with number of rows of B');
        end
        
        temp = p*eqn.D_ - eqn.K_;
        X1 = (p^2*eqn.M_ - temp)\(p*B(rowK+1:end,:) - B(1:rowK,:));
        X2 = (p*eqn.M_)\(B(1:rowK,:) - temp*X1);
        X  = [X1;X2];
        
    % implement solve (A+p*E)*X = B'
    case 'T'
        
        if (rowA ~= size(B,2))
            error('MESS:error_arguments','number of rows of A differs with number of columns of B');
        end
        
        temp = p*eqn.D_ -eqn.K_;
        X1 = (p^2*eqn.M_ - temp)\(p*B(:,rowK+1:end)' - B(:,1:rowK)');
        X2 = (p*eqn.M_)\(B(:,1:rowK)' - temp*X1);
        X  = [X1;X2];
    

end
