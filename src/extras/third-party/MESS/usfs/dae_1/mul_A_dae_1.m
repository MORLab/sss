function C = mul_A_dae_1(eqn, opts, opA, B, opB)

%% function mul_A perfoms operation C = opA(A_)*opB(B)
%
% Input:
%   eqn     structure contains field A_
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A_)
%           opA = 'N' performs A_*opB(B)
%           opA = 'T' performs A_'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opA(A_)*B
%           opB = 'T' performs opA(A_)*B'

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

%     A = [J1 J2;
%          J3 J4]   J4 regular
%
% Output:
% C = opA(A_)*opB(B)
%
%   uses size_dae_1

%% check input Paramters
if (~ischar(opA) || ~ischar(opB))
    error('MESS:error_arguments', 'opA or opB is not a char');
end

opA = upper(opA); opB = upper(opB);
if(~(opA == 'N' || opA == 'T'))
    error('MESS:error_arguments', 'opA is not ''N'' or ''T''');
end

if(~(opB == 'N' || opB == 'T'))
    error('MESS:error_arguments', 'opB is not ''N'' or ''T''');
end

if (~isnumeric(B)) || (~ismatrix(B))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(~isfield(eqn, 'A_')) || ~isnumeric(eqn.A_)
    error('MESS:error_arguments', 'field eqn.A_ is not defined');
end

if ~isfield(eqn, 'st')    || ~isnumeric(eqn.st)
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.')
end

n = size(eqn.A_,1);
st = eqn.st;
%% perfom multiplication
switch opA
    
    case 'N'
        switch opB
            
            %implement operation A_*B
            case 'N'
                if(st > size(B,1))
                    error('MESS:error_arguments', 'number of cols of A_ differs with rows of B');
                end
                C = eqn.A_(1 : st, 1 : st) * B(1 : st, : ) - eqn.A_(1 : st, st + 1 : n) * ...
                    (eqn.A_(st + 1 : n, st + 1 : n) \ (eqn.A_(st + 1 : n, 1 : st) * ...
                    B(1 : st, : )));
            
            %implement operation A_*B'
            case 'T'
                if(st > size(B, 2))
                    error('MESS:error_arguments', 'number of cols of A_ differs with cols of B');
                end
                C = eqn.A_(1 : st, 1 : st) * B( : , 1 : st)' - eqn.A_(1 : st, st + 1 : n) * ...
                    (eqn.A_(st + 1 : n, st + 1 : n) \ (eqn.A_(st + 1 : n, 1 : st) * ...
                    B( : , 1 : st)'));
        end
        
    case 'T'
        switch opB
            
            %implement operation A_'*B
            case 'N'
                if(st > size(B, 1))
                    error('MESS:error_arguments', 'number of rows of A_ differs with rows of B');
                end
                C = eqn.A_(1 : st, 1 : st)' * B(1 : st, : ) - eqn.A_(st + 1 : n, 1 : st)' * ...
                    (eqn.A_(st + 1 : n, st + 1 : n)' \ (eqn.A_(1 : st, st + 1 : n)' * ...
                    B(1 : st, : )));
                
            %implement operatio A_'*B'
            case 'T'
                if(st > size(B, 2))
                    error('MESS:error_arguments', 'number of rows of A_ differs with cols of B');
                end
                C = eqn.A_(1 : st, 1 : st)' * B( : , 1 : st)' - eqn.A_(st + 1 : n, 1 : st)' * ...
                    (eqn.A_(st + 1 : n, st + 1 : n)' \ (eqn.A_(1 : st, st + 1 : n)' * ...
                    B( : , 1 : st)'));
        end
        
end
end
