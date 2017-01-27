function X = sol_ApE_dae2_so(eqn, opts, opA, p, opE, B, opB)
%% function sol_ApE solves (opA(A_) + p*opE(E_))*X = opB(B) resp. performs X=(opA(A_)+p*opE(E_))\opB(B)
%
%
% A_ and E_ are assumed to be quadratic.
% Input:
%
%   eqn     structure contains A_ and E_
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A_)
%           opA = 'N' for A_
%           opA = 'T' for A_'
%
%   p       scalar Value
%
%   opE     character specifies the form of opE(E_)
%           opE = 'N' for E_
%           opE = 'T' for E_'
%
%   B       n-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' for B
%           opB = 'T' for B'
%
%   typeE   specifies whether E_ is Identity or not
%           typeE = 0 E_ is Identity
%           typeE = 1 E_ is not Identity
%
% Output
%
%   X       matrix fullfills equation (opA(A_)+p*opE(E_))*X = B
%

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

%% check input Paramters
if (~ischar(opA) || ~ischar(opE) || ~ischar(opB))
    error('MESS:error_arguments', 'opA, opE or opB is not a char');
end

opA = upper(opA); opE = upper(opE); opB = upper(opB);

if(~(opA == 'N' || opA == 'T'))
    error('MESS:error_arguments', 'opA is not ''N'' or ''T''');
end

if(~(opE == 'N' || opE == 'T'))
    error('MESS:error_arguments', 'opE is not ''N'' or ''T''');
end

if(~(opB == 'N' || opB == 'T'))
    error('MESS:error_arguments', 'opB is not ''N'' or ''T''');
end

if(~isnumeric(p))
   error('MESS:error_arguments','p is not numeric'); 
end

if (~isnumeric(B)) || (~ismatrix(B))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
for mat='MDKG'
    if(~isfield(eqn, sprintf('%c_',mat)) || ~eval(sprintf('isnumeric(eqn.%c_)',mat)))
        error('MESS:error_arguments', 'field eqn.%c_ is not defined',mat);
    end
end

nv = size(eqn.M_,1);
np = size(eqn.G_,1);

[rowB,colB] = size(B);

if(opB == 'N')
    if (rowB ~= 2*nv + np)
        B = [B; zeros(2 * nv + np - rowB, colB)];
    end
else
    if (colB ~= 2*nv + np)
        B = [B, zeros(rowB, 2 * nv + np - colB)];
    end
end

 switch opA
    
    case 'N'
      switch opE
        
        case 'N'
          
          switch opB
            %implement solve (A_+p*E_)*X=B
            case 'N'
              x23 = [eqn.K_-p*eqn.D_-p^2*eqn.M_, -p*eqn.G_';-p*eqn.G_,zeros(np,np)]\[eqn.K_*B(1:nv,:)-p*B(nv+1:2*nv,:);-p*B(2*nv+1:end,:)];
              X = [(B(1:nv,:)-x23(1:nv,:))./p;x23];
              %X = (A + p * E) \ B;
              
              %implement solve (A_+p*E_)*X=B'
            case 'T'
              x23 = [eqn.K_-p*eqn.D_-p^2*eqn.M_, -p*eqn.G_';-p*eqn.G_,zeros(np,np)]\[eqn.K_*B(:,1:nv)'-p*B(:,nv+1:2*nv)';-p*B(:,2*nv+1:end)'];
              X = [(B(:,1:nv)'-x23(1:nv,:))./p;x23];
              %X = (A + p * E) \ B';
          end
          
        case 'T'
          
          switch opB
            %implement solve (A_+p*E_)*X=B
            case 'N'
              x23 = [eqn.K_-p*eqn.D_-p^2*eqn.M_', -p*eqn.G_';-p*eqn.G_,zeros(np,np)]\[eqn.K_*B(1:nv,:)-p*B(nv+1:2*nv,:);-p*B(2*nv+1:end,:)];
              X = [(B(1:nv,:)-x23(1:nv,:))./p;x23];
              %X = (A + p * E) \ B;
              
              %implement solve (A_+p*E_)*X=B'
            case 'T'
              x23 = [eqn.K_-p*eqn.D_-p^2*eqn.M_', -p*eqn.G_';-p*eqn.G_,zeros(np,np)]\[eqn.K_*B(:,1:nv)'-p*B(:,nv+1:2*nv)';-p*B(:,2*nv+1:end)'];
              X = [(B(:,1:nv)'-x23(1:nv,:))./p;x23];
              %X = (A + p * E) \ B';
          end
                    
      end
      
    case 'T'
      switch opE
        
        case 'N'
             
            switch opB
                %implement solve (A_+p*E_)*X=B
                case 'N'
                  x23 = [eqn.K_'-p*eqn.D_'-p^2*eqn.M_, -p*eqn.G_';-p*eqn.G_,zeros(np,np)]\[B(1:nv,:)-p*B(nv+1:2*nv,:);-p*B(2*nv+1:end,:)];
                  X = [(B(1:nv,:)-eqn.K_'*x23(1:nv,:))./p;x23];
                  %X = (A + p * E) \ B;

                  %implement solve (A_+p*E_)*X=B'
                case 'T'
                  x23 = [eqn.K_'-p*eqn.D_'-p^2*eqn.M_, -p*eqn.G_';-p*eqn.G_,zeros(np,np)]\[B(:,1:nv)'-p*B(:,nv+1:2*nv)';-p*B(:,2*nv+1:end)'];
                  X = [(B(:,1:nv)'-eqn.K_'*x23(1:nv,:))./p;x23];
                  %X = (A + p * E) \ B';
            end
          
        case 'T'
          
            switch opB
                %implement solve (A_+p*E_)*X=B
                case 'N'
                  x23 = [eqn.K_'-p*eqn.D_'-p^2*eqn.M_', -p*eqn.G_';-p*eqn.G_,zeros(np,np)]\[B(1:nv,:)-p*B(nv+1:2*nv,:);-p*B(2*nv+1:end,:)];
                  X = [(B(1:nv,:)-eqn.K_'*x23(1:nv,:))./p;x23];
                  %implement solve (A_+p*E_)*X=B'
                case 'T'
                  x23 = [eqn.K_'-p*eqn.D_'-p^2*eqn.M_', -p*eqn.G_';-p*eqn.G_,zeros(np,np)]\[B(:,1:nv)'-p*B(:,nv+1:2*nv)';-p*B(:,2*nv+1:end)'];
                  X = [(B(:,1:nv)'-eqn.K_'*x23(1:nv,:))./p;x23];
                  %X = (A + p * E) \ B';
            end
      end
      
 end
 if opB == 'N'
     X = X(1 : rowB, :);
 else
     X = X(1 : colB, :);
 end
 % X = X(1 : 2*nv, :);
end

