function C = mul_A_dae_so_1(eqn, opts, opA, B, opB)

%% function mul_A perfoms operation C = opA(A_) * opB(B)
%
% Input:
%   eqn     structure contains field A_
%
%   opts    struct contains parameters for the algorithm
%s
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
%   uses no other dae_so_1 function

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
if (~isfield(eqn,'K_') || ~isnumeric(eqn.K_))
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.')
end
if (~isfield(eqn,'isSym'))
    isSym = 0;
else
    isSym = eqn.isSym;
end
if ~isfield(eqn, 'nd')    || ~isnumeric(eqn.nd)
    error('MESS:nd',...
    'Missing or Corrupted nd field detected in equation structure.');
end

nd = eqn.nd;

if(opB == 'N')
    rowB = size(B, 1);
else
    rowB = size(B, 2);
end

if(2 * nd ~= rowB)
    error('MESS:error_arguments', 'number of rows of B differs from number of cols of A ( 2 * nd)');
end

%% perfom multiplication
if isSym
    
    switch opB
        
        
        case 'N'
            C2 = eqn.K_(1 : nd, 1 : nd)* B(nd + 1 : end, :) ...
              - eqn.K_(1 : nd, nd + 1 : end) * ...
                (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                \ (eqn.K_(nd + 1 : end, 1 : nd) * B(nd + 1 : end, :)));
            C = [B(1 : nd, :); - C2];
            
            
        case 'T'
            C2 = eqn.K_(1 : nd, 1 : nd)* B( : , nd + 1 : end)'...
              - eqn.K_(1 : nd, nd + 1 : end) * ...
                ((eqn.K_(nd + 1 : end, nd + 1 : end) ...
                \ (eqn.K_(nd + 1 : end, 1 : nd) * B( : , nd + 1 : end)')));
            C = [B( : , 1 : nd)'; - C2];
    end
    
else
    switch opA
        
        case 'N'
            switch opB
                
                
                case 'N'
                    C2 = eqn.K_(1 : nd, 1 : nd) *B(nd + 1 : end, :) ...
                      - eqn.K_(1 : nd, nd + 1 : end) * ...
                        (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                        \ (eqn.K_(nd + 1 : end, 1 : nd) * B(nd + 1 : end, :)));
                    C = [B(1 : nd, :); - C2];
                    
                    
                case 'T'
                    C2 = eqn.K_(1 : nd, 1 : nd) * B( : , nd + 1 : end)'...
                      - eqn.K_(1 : nd, nd + 1 : end) * ...
                        (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                        \ (eqn.K_(nd + 1 : end, 1 : nd) * B( : , nd + 1 : end)'));
                    C = [B( : , 1 : nd)'; - C2];
            end
            
        case 'T'
            switch opB
                
                
                case 'N'
                    C2 = eqn.K_(1 : nd, 1 : nd)'* B(nd + 1 : end, :) ...
                      - eqn.K_(nd + 1 : end, 1 : nd)' * ...
                        (eqn.K_(nd + 1 : end, nd + 1 : end)' ...
                        \ (eqn.K_(1 : nd, nd + 1 : end)' * B(nd + 1 : end, :)));
                    C = [B(1 : nd, :); - C2];
                    
                    
                case 'T'
                    C2 = eqn.K_(1 : nd, 1 : nd)'* B( : , nd + 1 : end)'...
                      - eqn.K_(nd + 1 : end, 1 : nd)' * ...
                        (eqn.K_(nd + 1 : end, nd + 1 : end)' ...
                        \ (eqn.K_(1 : nd, nd + 1 : end)' * B( : , nd + 1 : end)'));
                    C = [B( : , 1 : nd)'; - C2];
            end
            
    end
end
end
