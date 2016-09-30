function X=sol_ApE_solveLse(eqn, opts,opA,p,opE,C,opC)

% function X=sol_ApE_default(eqn, opts,opA,p,opE,C,opC)
%
% This function returns X = (A_ + p*E_)\C, where matrices A_ and E_ given by structure eqn and input matrix C could be transposed.
% Matrices A_ and E_ are assumed to be quadratic.
%
%   Inputs:
%
%   eqn     structure containing fields 'A_' and 'E_'
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' solves (A_ + p* opE(E_))*X = opC(C) 
%           opA = 'T' solves (A_' + p* opE(E_))*X = opC(C) 
%   p       scalar value
%   opE     character specifying the shape of E_
%           opE = 'N' solves (opA(A_) + p* E_)*X = opC(C) 
%           opE = 'T' solves (opA(A_) + p* E_')*X = opC(C) 
%   C       n-x-p matrix
%   opC     character specifies the form of opC(C)
%           opC = 'N' solves (opA(A_) + p* opE(E_))*X = C
%           opC = 'T' solves (opA(A_) + p* opE(E_))*X = C'
%
%   Output:
%
%   X       matrix fullfilling equation (opA(A_)+p*opE(E_))*X = opC(C)
%
% This function does not use other default functions.

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
% Copyright (C) Jens Saak, Martin Koehler and others 
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

if(~isnumeric(p)) || (length(p) ~= 1)
    error('MESS:error_arguments','p is not numeric');
end

if (~isnumeric(C)) || (~ismatrix(C))
    error('MESS:error_arguments','C has to ba a matrix');
end

%% check data in eqn structure
if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
if(eqn.haveE ==1)
    if(~isfield(eqn,'E_') || ~isfield(eqn,'A_'))
        error('MESS:error_arguments','field eqn.E_ or eqn.A_ is not defined');
    end
else
    if(~isfield(eqn,'A_'))
        error('MESS:error_arguments','field eqn.A_ is not defined');
    end    
end

[rowA,colA] = size(eqn.A_);

if real(p)>=0
  error('Real parts of entries of p must be negative!');
end

%% check if new LU is necessary
% get information about last solveLse call
lastLse=last_solveLse();

% check if newlu is necessary
if ~isfield(lastLse,'p') || isempty(lastLse.p) || ~isfield(lastLse,'opA') || isempty(lastLse.opA)...
       || ~isfield(lastLse,'opE') || isempty(lastLse.opE) ||  ~isfield(lastLse,'solveLse') || isempty(lastLse.solveLse)
    newlu=true;
elseif p==lastLse.p && opA==lastLse.opA && opE == lastLse.opE && strcmp(lastLse.solveLse,'sol_ApE')
    newlu=false;
else
    newlu=true;
end

if newlu
    jCol=1; % new lu
    s0=ones(1,size(C,2)+1)*p; % keep LU factors persistent in solveLse (jCol<length(s0))
    Rt=[speye(size(C,2)),[1;zeros(size(C,2)-1,1)]];
else
    jCol=2; % reuse old lu
    s0=ones(1,size(C,2)+2)*p; % keep LU factors persistent in solveLse (jCol<length(s0))
    Rt=[[zeros(size(C,2)-1,1);1],speye(size(C,2)),[1;zeros(size(C,2)-1,1)]];
end

lastLse.p=p;
lastLse.opA=opA;
lastLse.opB=[];
lastLse.opC=opC;
lastLse.opE=[];
lastLse.solveLse='sol_ApE';

% update lastLse
last_solveLse(lastLse);

  
switch opA

  case 'N'

      switch opC

          %implement solve (A_+p*I)*X=C
          case 'N'

              if(rowA~=size(C,1)) 
                  error('MESS:error_arguments','number of rows of A differs with number of rows of C');
              end
              X=zeros(size(C,1),size(s0,2));
              for i=1:size(C,2)
                X = solveLse(jCol, X, eqn.A_, C, eqn.E_, -s0, Rt, opts.solveLse);
                jCol=jCol+1;
              end
              if length(s0)==size(C,2)+1
                  X=X(:,1:size(C,2));
              else
                  X=X(:,2:size(C,2)+1);
              end        

          %implement solve (A_+p*I)*X=C'
          case 'T'

              if(rowA~=size(C,2)) 
                  error('MESS:error_arguments','number of rows of A differs with number of columns of C');
              end
              X=zeros(size(C',1),size(s0,2));
              for i=1:size(C,2)
                  X = solveLse(jCol, X, eqn.A_, C', eqn.E_, -s0, Rt, opts.solveLse);
                  jCol=jCol+1;
              end
              if length(s0)==size(C,2)+1
                  X=X(:,1:size(C,2));
              else
                  X=X(:,2:size(C,2)+1);
              end     

      end

  case 'T'

      switch opC

          %implement solve (A_'+p*I)*X=C
          case 'N'

              if(colA~=size(C,1)) 
                  error('MESS:error_arguments','number of columns of A differs with number of rows of C');
              end
              X=zeros(size(C,1),size(s0,2));
              for i=1:size(C,2)
                  X = solveLse(jCol, X, eqn.A_', C, eqn.E_, -s0, Rt, opts.solveLse);
                  jCol=jCol+1;
              end
              if length(s0)==size(C,2)+1
                  X=X(:,1:size(C,2));
              else
                  X=X(:,2:size(C,2)+1);
              end     

          %implement solve (A_'+p*I)*X=C'
          case 'T'

              if(colA~=size(C,2)) 
                  error('MESS:error_arguments','number of columns of A differs with number of columns of C');
              end
              X=zeros(size(C',1),size(s0,2));
              for i=1:size(C,2)
                  X = solveLse(jCol, X, eqn.A_', C', eqn.E_, -s0, Rt, opts.solveLse);
                  jCol=jCol+1;
              end
              if length(s0)==size(C,2)+1
                  X=X(:,1:size(C,2));
              else
                  X=X(:,2:size(C,2)+1);
              end     

      end

end


end
