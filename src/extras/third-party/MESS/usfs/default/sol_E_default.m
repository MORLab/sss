function X=sol_E_default(eqn, opts,opE,B,opB)

% function X=sol_E_default(eqn, opts,opE,B,opB)
%
% This function returns X = E_\B, where matrix E_ given by structure eqn and input matrix B could be transposed.
% Matrix E_ is assumed to be quadratic and has the same size as A_ in structure eqn.
%
%   Inputs:
%
%   eqn       structure containing field 'E_'
%   opts      structure containing parameters for the algorithm
%   opE       character specifying the shape of E_
%                  opE = 'N' solves E_*X = opB(B)
%                  opE = 'T' solves E_'*X = opB(B)  
%   B         p-x-q matrix
%   opB       character specifying the shape of B
%                  opB = 'N' solves  opE(E_)*X = B     
%                  opB = 'T' solves  opE(E_)*X = B' 
%
%   Output:
%
%   X         matrix fullfilling equation  opE(E_)*X = opB(B)
%
% This function uses another default function size_default(eqn, opts) to obtain the number of rows of matrix A_ in structure eqn,
% that should be equal to the number of rows of matrix E_.
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
if(~isfield(eqn,'E_'))
    error('MESS:error_arguments','field eqn.E_ is not defined');
end

rowE = size_default(eqn, opts);
colE = rowE; %we only consider square systems

%% perform solve operations
switch opE
    
    case 'N'
        switch opB
            
            %implement solve E_*X=B
            case 'N'
                if(rowE~=size(B,1))
                    error('MESS:error_arguments','number of rows of E_ differs with number of rows of B');
                end
                X = eqn.E_\B;
            
            %implement solve E_*X=B'
            case 'T'
                if(rowE~=size(B,2))
                    error('MESS:error_arguments','number of rows of E_ differs with number of columns of B');
                end
                X = eqn.E_\B';
        end
        
    case 'T'
        switch opB
            
            %implement solve E_'*X=B
            case 'N'
                if(colE~=size(B,1))
                    error('MESS:error_arguments','number of columns of E_ differs with number of rows of B');
                end
                X = eqn.E_'\B;
                
            %implement solve E_'*X=B'
            case 'T'
                if(colE~=size(B,2))
                    error('MESS:error_arguments','number of columns of E_ differs with number of columns of B');
                end
                X = eqn.E_'\B';
        end
        
end

end

