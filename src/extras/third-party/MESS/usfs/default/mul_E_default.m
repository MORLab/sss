function C=mul_E_default(eqn, opts,opE,B,opB)

% function C=mul_E_default(eqn, opts,opE,B,opB)
%
% This function returns C = E_*B, where matrix E_ given by struture eqn and input matrix B could be transposed.
% Matrix E_ is assumed to be quadratic and has the same size as A_ in structure eqn.
%
%   Inputs:
%
%   eqn     structure containing field 'E_'
%   opts    structure containing parameters for the algorithm
%   opE     character specifying the shape of E_
%           opE = 'N' performs E_*opB(B)
%           opE = 'T' performs E_'*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs opE(E_)*B
%           opB = 'T' performs opE(E_)*B'
%
%   Output:
%
%   C = opE(E_)*opB(B)
%
% This function uses another default function size_default(eqn,opts) to obtain the number of rows of matrix A_ in structure eqn,
% that should be equal to the number of rows of the matrix E_.

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
colE =rowE;

%% perform multiplication
switch opE
    
    case 'N'
        switch opB
            
            %implement operation E_*B
            case 'N'
                if(colE~=size(B,1))
                    error('MESS:error_arguments','number of columns of E_ differs with number of rows of B');
                end
                C = eqn.E_*B;
            
            %implement operation E_*B'
            case 'T'
                if(colE~=size(B,2))
                    error('MESS:error_arguments','number of columns of E_ differs with number of columns of B');
                end
                C = eqn.E_*B';
        end
        
    case 'T'
        switch opB
            
            %implement operation E_'*B
            case 'N'
                if(rowE~=size(B,1))
                    error('MESS:error_arguments','number of rows of E_ differs with number of rows of B');
                end
                C = eqn.E_'*B;
                
            %implement operation E_'*B'
            case 'T'
                if(rowE~=size(B,2))
                    error('MESS:error_arguments','number of rows of E_ differs with number of columns of B');
                end
                C = eqn.E_'*B';
        end
        
end

end
