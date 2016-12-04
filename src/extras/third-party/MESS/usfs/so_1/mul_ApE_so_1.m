function C=mul_ApE_so_1(eqn, opts,opA,p,opE,B,opB)

% function C=mul_ApE_so_1(eqn, opts,opA,p,opE,B,opB)
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
% This function returns C = (A+p*E)*B, where matrices A and E given by structure eqn and input matrix B could be transposed.
% Matrix A is assumed to be quadratic and has a size of 2*size(K).
%
%   Inputs:
%
%   eqn     structure containing data for matrices A (fields 'K_' and 'D_') and E (fields 'K_' and 'M_')
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
%       (| 0 -K|      |-K  0|)
%   C = (|-K -D| + p* | 0  M|)*opB(B)
%
% This function does not use other so1 functions.
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
    if~isfield(eqn,'D_') || ~isnumeric(eqn.D_) ||...
            ~isfield(eqn,'K_') || ~isnumeric(eqn.K_)
        error('MESS:error_arguments',...
            'field eqn.K_ or eqn.D_ is not defined or corrupted');
    end
end

[rowK, colK] = size(eqn.K_);
rowA = 2*rowK;
colA = 2*colK;



if(eqn.haveE==1)
%% perform operations for E ~= Identity
    switch opB
        
        %implement multiplication (A+p*E)*B=C
        case 'N'
            if(colA~=size(B,1))
                error('MESS:error_arguments','number of columns of A differs with number of rows of B');
            end
            C = [(-p*(eqn.K_*B(1:rowK,:)) -(eqn.K_*B(rowK+1:end,:)));
                 (-(eqn.K_*B(1:rowK,:)) + p*(eqn.M_*B(rowK+1:end,:)) - eqn.D_*B(rowK+1:end,:))];
            
        %implement multiplication (A+p*E)*B'=C
        case 'T'
            if(colA~=size(B,2))
                error('MESS:error_arguments','number of columns of A differs with number of columns of B');
            end
            C = [(-p*(eqn.K_*B(:,1:colK)') -(eqn.K_*B(:,colK+1:end)'));...
                 (-(eqn.K_*B(:,1:colK)')   +  p*(eqn.M_*B(:,colK+1:end)') - eqn.D_*B(:,colK+1:end)')];
    end
    
    
elseif(eqn.haveE==0)
%% perform operations for E = Identity
    
    switch opB
    
        %implement multiplication (A+p*I)*B=C
        case 'N'
            if(colA~=size(B,1))
                error('MESS:error_arguments','number of columns of A differs with number of rows of B');
            end
            C = [(p*B(1:rowK,:)      -(eqn.K_*B(rowK+1:end,:)));...
                 (-(eqn.K_*B(1:rowK,:))+  p*B(rowK+1:end,:) - eqn.D_*B(rowK+1:end,:))];
            
        %implement multiplication (A+p*I)*B'=C
        case 'T'
            if(colA~=size(B,2))
                error('MESS:error_arguments','number of columns of A differs with number of columns of B');
            end
            C = [(p*B(:,1:colK)'      -(eqn.K_*B(:,colK+1:end)'));...
                 (-(eqn.K_*B(:,1:colK)')+  p*B(:,colK+1:end)' - eqn.D_*B(:,colK+1:end)')];
    end
    
end
end
