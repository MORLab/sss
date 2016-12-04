function C = mul_E_dae_2(eqn, opts, opE, B, opB)

%% function mul_A perfoms operation C = opE(E_)*opB(B)
%
% Input:
%   eqn     structure contains field E_
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E_)
%           opE = 'N' performs E_*opB(B)
%           opE = 'T' performs E_'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opE(E_)*B
%           opB = 'T' performs opE(E_)*B'
%
% Output:
% C = opE(E_)*opB(B)
%
%   uses no other dae_1 function

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

if(~isfield(eqn, 'S_')) || ~isnumeric(eqn.S_)
    error('MESS:error_arguments', ['field eqn.S_ is not defined. Did ' ...
                        'you forget to run mul_E_pre?']);
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

%% perfom multiplication
switch opE
    
    case 'N'
        switch opB
            
            %implement operation E_*B
            case 'N'
              C = eqn.S_(1 : rowB, 1 : rowB) * B;
            
            %implement operation E_*B'
            case 'T'
              C = eqn.S_(1 : rowB, 1 : rowB) * B';
        end
        
    case 'T'
        switch opB
            
            %implement operation E_'*B
            case 'N'
                C = eqn.S_(1 : rowB, 1 : rowB)' * B;
                
            %implement operation E_'*B'
            case 'T'
                C = eqn.S_(1 : rowB, 1 : rowB)' * B';
        end
        
end

end
