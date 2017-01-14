function X = sol_A_dae3_so(eqn, opts, opA, B, opB)
%  function sol_A solves solves opA(A_)*X = opB(B)
%
% Input:
%  eqn       structure with field A_
%  opts      struct contains parameters for the algorithm
%  opA       character specifies the form of opA(A_) 
%                  opA = 'N' solves A_*X=opB(B)
%                  opA = 'T' solves A_^T*X=opB(B)
%  B         p-x-q matrix
%  opB       character specifies the form of opB(B)
%                  opB = 'N' solves A_*X=B
%                  opB = 'T' solves A_*X=B^T 
%
% Output:
%  X       matrix fullfills equation opA(A_)X = opB(B)
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
if (~ischar(opA) || ~ischar(opB))
    error('MESS:error_arguments', 'opA or opB is not a char');
end

opA = upper(opA); opB = upper(opB);
if(~(opA == 'N' || opA == 'T'))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(~(opB == 'N' || opB == 'T'))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if (~isnumeric(B)) || (~ismatrix(B))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
for mat='DKG'
    if(~isfield(eqn, sprintf('%c_',mat)) || ~eval(sprintf('isnumeric(eqn.%c_)',mat)))
        error('MESS:error_arguments', 'field eqn.%c_ is not defined',mat);
    end
end

nv = size(eqn.M_,1);
np = size(eqn.G_,1);

%% solve

if (opB=='N' && (size(B,1)==(2*nv+np))) || (opB=='T' && (size(B,2)==(2*nv+np)))
switch opA
    
    case 'N'
        switch opB
            
            %implement solve A_*X=B
            case 'N'
                x = [eqn.K_, eqn.G_';eqn.G_,sparse(np,np)]\ [B(nv+1:2*nv,:)-eqn.D_*B(1:nv,:);B(2*nv+1:end,:)];
                X = [x(1:nv,:);B(1:nv,:);x(nv+1:end,:)];
            
            %implement solve A_*X=B'
            case 'T'
                x = [eqn.K_, eqn.G_';eqn.G_,sparse(np,np)]\ [B(:,nv+1:2*nv)'-eqn.D_*B(:,1:nv)';B(:,2*nv+1:end)'];
                X = [x(1:nv,:);B(:,1:nv)';x(nv+1:end,:)];
        end
        
    case 'T'
        switch opB
             
            %implement solve A_'*X=B
            case 'N'
                x = [eqn.K_',eqn.G_';eqn.G_,sparse(np,np)]\[B(1:nv,:);B(2*nv+1:end,:)];
                X = [B(nv+1:2*nv,:)-eqn.D_'*x(1:nv,:);x(1:nv,:);x(nv+1:end,:)];
                
            %implement solve A_'*X=B'
            case 'T'
                x = [eqn.K_',eqn.G_';eqn.G_,sparse(np,np)]\[B(:,1:nv)';B(:,2*nv+1:end)'];
                X = [B(:,nv+1:2*nv)'-eqn.D_'*x(1:nv,:);x(1:nv,:);x(nv+1:end,:)];
        end
        
end
elseif (opB=='N' && (size(B,1)==(2*nv))) || (opB=='T' && (size(B,2)==(2*nv)))
    error('MESS:error_usage','mul_A_dae2_so is only coded for shift parameter computation');
else
    error('MESS:error_arguemnts', 'B has wrong number of cols');
end



end
