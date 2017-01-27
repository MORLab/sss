function X = sol_E_dae_2(eqn, opts, opE, B, opB) 
%% function sol_E_dae_so_1 solves opE(E)*X = opB(B) resp. performs X=opE(E)\opB(B)
%
% Input:
%   eqn     structure contains data for E (M_,K_)
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' solves E *X = opB(B)
%           opE = 'T' sovles E'*X = opB(B)
%
%   B       p-x-q matrix 
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' solves opE(E)*X = B
%           opB = 'T' solves opE(E)*X = B'
%
% Output
%
%   X       matrix fullfills equation opE(E)*X = opB(B)
%
%   uses no other dae_1 function
%% check input Paramters

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
if (~ischar(opE) || ~ischar(opB))
    error('MESS:error_arguments', 'opE or opB is not a char');
end

opE = upper(opE); opB = upper(opB);
if(~(opE == 'N' || opE == 'T'))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(~(opB == 'N' || opB == 'T'))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end
if (~isnumeric(B)) || (~ismatrix(B))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(~isfield(eqn, 'E_')) || ~isnumeric(eqn.E_)
    error('MESS:error_arguments', 'field eqn.E_ is not defined');
end
if(~isfield(eqn, 'S_'))
    error('MESS:error_arguments', ['field eqn.S_ is not defined. Did ' ...
                        'you forget to run sol_E_pre?']);
end
if ~isfield(eqn, 'st')    || ~isnumeric(eqn.st)
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.')
end

st = eqn.st;
switch opB
  case 'N'
    rowB=size(B,1);
  case 'T'
    rowB=size(B,2);    
end
if rowB~=st && rowB~=size(eqn.E_,1) 
  error('MESS:error_arguments', 'size of B does not match data in E');
end


%% solve
switch opE
    
    case 'N'
        switch opB
            
            %implement solve E*X=B
            case 'N'
                if( size(eqn.S_,1) ~= size(B, 1))
                    error('MESS:error_arguments','number of rows of E differs with rows of B');
                end
                X = eqn.S_ \ B;
            
            %implement solve A*X=B'
            case 'T'
                if( size(eqn.S_,1) ~= size(B, 2))
                    error('MESS:error_arguments','number of rows of E differs with cols of B');
                end
                X = eqn.S_ \ B';
        end
        
    case 'T'
        switch opB
            
            %implement solve E'*X=B
            case 'N'
                if( size(eqn.S_,1) ~= size(B, 1))
                    error('MESS:error_arguments','number of cols of A_ differs with rows of B');
                end
                X = eqn.S_' \ B;
                
            %implement solve A_'*X=B'
            case 'T'
                if( size(eqn.S_,1) ~= size(B, 2))
                    error('MESS:error_arguments','number of cols of A_ differs with cols of B');
                end
                X = eqn.S_' \ B';
        end
        
end

end
