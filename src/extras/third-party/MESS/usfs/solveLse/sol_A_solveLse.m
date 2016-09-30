function X=sol_A_solveLse(eqn, opts,opA,B,opB)

% function X=sol_A_default(eqn, opts,opA,B,opB)
%
% This function returns X = A_\B, where matrix A_ given by structure eqn and input matrix B could be transposed.
% Matrix A_ is assumed to be quadratic.
% 
%    Inputs:
%
%    eqn       structure containing field 'A_'
%    opts      structure containing parameters for the algorithm
%    opA       character specifying the shape of A_
%                  opA = 'N' solves  A_*X = opB(B) 
%                  opA = 'T' solves  A_'*X = opB(B) 
%    B         p-x-q matrix
%    opB       character specifying the shape of B
%                  opB = 'N' solves  opA(A_)*X = B    
%                  opB = 'T' solves  opA(A_)*X = B'
%
%    Output:
%
%    X         matrix fullfilling equation  opA(A_)*X = opB(B) 
%
% This function uses another default function size_default(eqn, opts) to obtain the number of rows of matrix A_ in structure eqn.
%
% Maximilian Behr 2012/11/1

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
if (~ischar(opA) || ~ischar(opB))
    error('MESS:error_arguments', 'opA or opB is not a char');
end

opA = upper(opA); opB = upper(opB);
if(~(opA=='N' || opA=='T'))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(~(opB=='N' || opB=='T'))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if (~isnumeric(B)) || (~ismatrix(B))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(~isfield(eqn,'A_'))
    error('MESS:error_arguments','field eqn.A_ is not defined');
end

rowA = size_default(eqn, opts);
colA = rowA;

%% check if new LU is necessary

% get information about last solveLse call
lastLse=last_solveLse();

% check if newlu is necessary
if ~isfield(lastLse,'opA') || isempty(lastLse.opA) || ~isfield(lastLse,'opB')...
    || isempty(lastLse.opB) ||  ~isfield(lastLse,'solveLse') || isempty(lastLse.solveLse)
    newlu=true;
elseif opA==lastLse.opA && opB == lastLse.opB && strcmp(lastLse.solveLse,'sol_A')
    newlu=false;
else
    newlu=true;
end

if newlu
    jCol=1; % new lu
    s0=ones(1,size(B,2)+1)*0; % keep LU factors persistent in solveLse (jCol<length(s0))
    Rt=[speye(size(B,2)),[1;zeros(size(B,2)-1,1)]];
else
    jCol=2; % reuse old lu
    s0=ones(1,size(B,2)+2)*0; % keep LU factors persistent in solveLse (jCol<length(s0))
    Rt=[[zeros(size(B,2)-1,1);1],speye(size(B,2)),[1;zeros(size(B,2)-1,1)]];
end



lastLse.p=0;
lastLse.opA=opA;
lastLse.opB=opB;
lastLse.opC=[];
lastLse.opE=[];
lastLse.solveLse='sol_A';

% update lastLse
last_solveLse(lastLse);


%% perform solve operations
switch opA
    
    case 'N'
        switch opB
            
            %implement solve A_*X=B
            case 'N'
                if(rowA~=size(B,1))
                    error('MESS:error_arguments','number of rows of A_ differs with number of rows of B');
                end
                X=zeros(size(B,1),size(s0,2));
                for i=1:size(B,2)
                    X = solveLse(jCol, X, eqn.A_,B, eqn.E_, s0, Rt, opts.solveLse);
                    jCol=jCol+1;
                end
                if length(s0)==size(B,2)+1
                    X=X(:,1:size(B,2));
                else
                    X=X(:,2:size(B,2)+1);
                end 
            
            %implement solve A_*X=B'
            case 'T'
                if(rowA~=size(B,2))
                    error('MESS:error_arguments','number of rows of A_ differs with number of columns of B');
                end
                X=zeros(size(B',1),size(s0,2));
                for i=1:size(B,2)
                    X = solveLse(jCol, X, eqn.A_,B', eqn.E_, s0, Rt, opts.solveLse);
                    jCol=jCol+1;
                end
                if length(s0)==size(B,2)+1
                    X=X(:,1:size(B,2));
                else
                    X=X(:,2:size(B,2)+1);
                end 
        end
        
    case 'T'
        switch opB
            
            %implement solve A_'*X=B
            case 'N'
                if(colA~=size(B,1))
                    error('MESS:error_arguments','number of columns of A_ differs with number of rows of B');
                end
                X=zeros(size(B,1),size(s0,2));
                for i=1:size(B,2)
                    X = solveLse(jCol, X, eqn.A_',B, eqn.E_, s0, Rt, opts.solveLse);
                    jCol=jCol+1;
                end
                if length(s0)==size(B,2)+1
                    X=X(:,1:size(B,2));
                else
                    X=X(:,2:size(B,2)+1);
                end 
                
            %implement solve A_'*X=B'
            case 'T'
                if(colA~=size(B,2))
                    error('MESS:error_arguments','number of columns of A_ differs with number of columns of B');
                end
                X=zeros(size(B',1),size(s0,2));
                for i=1:size(B,2)
                    X = solveLse(jCol, X, eqn.A_',B', eqn.E_, s0, Rt, opts.solveLse);
                    jCol=jCol+1;
                end
                if length(s0)==size(B,2)+1
                    X=X(:,1:size(B,2));
                else
                    X=X(:,2:size(B,2)+1);
                end 
                
        end
        
end



end

