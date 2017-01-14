function C=mul_ApE_so_2(eqn, opts,opA,p,opE,B,opB)

% function C=mul_ApE_so_2(eqn, opts,opA,p,opE,B,opB)
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
% This function returns C = (A+p*E)*B, where matrices A and E given by structure eqn and input matrix B could be transposed.
% Matrix A is assumed to be quadratic and has a size of 2*size(K).
%
%   Inputs:
%
%   eqn     structure containing data for matrices A (fields 'K_' and 'M_') and E (fields 'D_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' performs (A + p*opE(E))*opB(B)
%           opA = 'T' performs (A' + p*opE(E))*opB(B)
%   p       scalar value
%   opE     character specifying the shape of E
%           opA = 'N' performs (opA(A) + p*E)*opB(B)
%           opA = 'T' performs (opA(A) + p*E')*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs (opA(A) + p*opE(E))*B
%           opB = 'T' performs (opA(A) + p*opE(E))*B'
%
%   Output:
%
%       (|-K  0|      | D  M|)
%   C = (| 0  M| + p* | M  0|)*opB(B)
%
% This function does not use other so3 functions.
%
% ATTENTION: opA,opE are not used since matrices A and E are symmetric. 

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
if (~ischar(opA) || ~ischar(opB) || ~ischar(opE))
    error('MESS:error_arguments', 'opA, opB or opE is not a char');
end

opA = upper(opA); opB = upper(opB); opE = upper(opE);

if(~(opA=='N' || opA=='T'))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(~(opB=='N' || opB=='T'))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if(~(opE=='N' || opE=='T'))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
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

[rowK, colK] = size(eqn.K_);
rowA = 2*rowK;
colA = 2*colK;

%% perform multiplication
switch opB
    

    % implement operation (A+p*E)*B = C
    case 'N' 
        
        if(colA ~= size(B,1))
            error('MESS:error_arguments','number of columns of A differs with number of rows of B');
        end
        temp = p*eqn.M_;
        C = [ (p*eqn.D_ -eqn.K_)*B(1:rowK,:) + temp*B(rowK+1:end,:); ...
              temp*B(1:rowK,:) + eqn.M_*B(rowK+1:end,:)];
        
    % implement operation (A+p*E)*B'= C
    case 'T'
        
        if(colA ~= size(B,2))
            error('MESS:error_arguments','number of columns of A differs with number of columns of B');
        end
        
        temp = p*eqn.M_;
        C = [(p*eqn.D_ - eqn.K_)*B(:,1:rowK)' + temp*B(:,rowK+1:end)';...
             temp*B(:,1:rowK)' + eqn.M_*B(:,rowK+1:end)'];

end

end
