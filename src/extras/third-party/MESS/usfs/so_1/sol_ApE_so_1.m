function X=sol_ApE_so_1(eqn, opts,opA,p,opE,C,opC)

% function X=sol_ApE_so_1(eqn, opts,opA,p,opE,C,opC)
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
% This function returns X =(A+p*E)\C, where matrices A and E given by structure eqn and input matrix C could be transposed.
% Matrices A and E are assumed to be quadratic.
%
%   Inputs:
%
%   eqn     structure containing data for matrices A (fields 'K_' and 'D_') and E (fields 'K_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' solves (A + p* opE(E))*X = opC(C) 
%           opA = 'T' solves (A' + p* opE(E))*X = opC(C) 
%   p       scalar value
%   opE     character specifying the shape of E
%           opE = 'N' solves (opA(A) + p* E)*X = opC(C) 
%           opE = 'T' solves (opA(A) + p* E')*X = opC(C) 
%   C       n-x-p matrix
%   opC     character specifies the shape of C
%           opC = 'N' solves (opA(A) + p* opE(E))*X = C
%           opC = 'T' solves (opA(A) + p* opE(E))*X = C'
%
%   Output:
%
%                                       (| 0 -K|     |-K  0|)
%   X       matrix fullfilling equation (|-K  D| + p*| 0  M|)*X = opC(C)
%
%   This function does not use other so1 functions.
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
if (~ischar(opA) || ~ischar(opE) || ~ischar(opC))
    error('MESS:error_arguments', 'opA, opE or opC is not a char');
end

opA = upper(opA); opE = upper(opE); opC = upper(opC);

if(~(opA=='N' || opA=='T'))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(~(opE=='N' || opE=='T'))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(~(opC=='N' || opC=='T'))
    error('MESS:error_arguments','opC is not ''N'' or ''T''');
end

if(~isnumeric(p))
    error('MESS:error_arguments','p is not numeric');
end

if (~isnumeric(C)) || (~ismatrix(C))
    error('MESS:error_arguments','C has to ba a matrix');
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
%% perform solve operations for E ~= Identity   

    switch opC

        %implement solve (A+p*E)*X=C
        case 'N'
            
            if (rowA ~= size(C,1))
                error('MESS:error_arguments','number of rows of A differs with number of rows of C');
            end
            
            X2=(p^2*eqn.M_-eqn.D_*p+eqn.K_)\(p*C(rowK+1:end,:)-C(1:rowK,:));
            X1=(-p*eqn.K_)\(C(1:rowK,:)+eqn.K_*X2);
            X=[X1;X2];
            
        %implement solve (A+p*E)*X=C'
        case 'T'
            
            if (rowA ~= size(C,2))
                error('MESS:error_arguments','number of rows of A differs with number of columns of C');
            end
            
            X2=(p^2*eqn.M_-eqn.D_*p+eqn.K_)\(p*C(:,rowK+1:end)'-C(:,1:rowK)');
            X1=(-p*eqn.K_)\(C(:,1:rowK)'+eqn.K_*X2);
            X=[X1;X2];
    end
    
elseif(eqn.haveE==0)
%% perform solve operations for E = Identity
    
    switch opC
        
                %implement solve (A+p*I)*X=C
                case 'N'
                    
                    if (rowA ~= size(C,1))
                        error('MESS:error_arguments','number of rows of A differs with number of rows of C');
                    end
                    
                    X2=(-(eqn.K_^2) -p*eqn.D_ +p^2*speye(rowK,colK))\(p*C(rowK+1:end,:)+eqn.K_*C(1:rowK,:));
                    X1=(C(1:rowK,:)+eqn.K_*X2)/p;                    
                    X=[X1;X2];
        
                %implement solve (A+p*I)*X=C'
                case 'T'
                    
                    if (rowA ~= size(C,2))
                        error('MESS:error_arguments','number of rows of A differs with number of columns of C');
                    end
                    
                    X2=(-(eqn.K_^2) -p*eqn.D_ +p^2*speye(rowK,colK))\(p*C(:,rowK+1:end)'+eqn.K_*C(:,1:rowK)');
                    X1=(C(:,1:rowK)'+eqn.K_*X2)/p;                    
                    X=[X1;X2];      
    end
    
end
end

