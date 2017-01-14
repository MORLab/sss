function C=mul_E_so_1(eqn, opts,opE,B,opB)

% function C=mul_E_so_1(eqn, opts,opE,B,opB)
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
% This function returns C = E*B, where matrix E given by structure eqn and input matrix B could be transposed.
% Matrix E is assumed to be quadratic and has a same size of 2*size(K).
%
%   Inputs:
%
%   eqn     structure containing data for matrix E (fields 'K_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opE     character specifying the shape of E
%           opE = 'N' performs E*opB(B)
%           opE = 'T' performs E'*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs opE(E)*B
%           opB = 'T' performs opE(E)*B'
%
%   Output:
%
%       |-K  0|
%   C = | 0  M| *opB(B)
%
% This function does not use other so1 functions. 
%
% ATTENTION: opE is not used since matrix E is symmetric

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
if(~isfield(eqn,'K_') || ~isnumeric(eqn.K_) || ~isfield(eqn,'M_') ...
        || ~isnumeric(eqn.M_))
    error('MESS:error_arguments',...
        'E consists of K and M, field eqn.K_ or eqn.M_ is not defined');
end


[rowK, colK] = size(eqn.K_);
rowE = 2*rowK;
colE = 2*colK;


%% perform multiplication
switch opB
    
    %implement operation E*B
    case 'N'
        if(colE~=size(B,1))
            error('MESS:error_arguments','number of columns of E differs with number of rows of B');
        end
        C = [-eqn.K_*B(1:rowK,:);
              eqn.M_*B(rowK+1:end,:)];
        
    %implement operation E*B'
    case 'T'
        if(colE~=size(B,2))
            error('MESS:error_arguments','number of columns of E differs with number of columns of B');
        end
        C=[-eqn.K_*B(:,1:colK)';...
            eqn.M_*B(:,colK+1:end)'];                    
         
end


end

