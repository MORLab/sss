function X = sol_A_dae_so_1(eqn, opts, opA, B, opB)

%% function sol_A solves opA(A) * X = opC(B) resp. performs X = opA(A) \ opB(B)
%
%
% M, D, K are assumed to be quadratic.
% Input:
%
%   eqn     structure contains  data for A (D_,K_) and E (M_,K_)
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A)
%           opA = 'N' for A
%           opA = 'T' for A'
%
%   B       n-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' for B
%           opB = 'T' for B'
%
% Output
%
%   X       matrix fullfills equation opA(A) * X = B
%
%   uses no other dae_so_1 function

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

rowK = size(eqn.K_, 1);
nd = eqn.nd;
na = rowK - nd;

if(opB == 'N')
    rowB = size(B, 1);
    colB = size(B, 2);
else
    rowB = size(B, 2);
    colB = size(B, 1);
end

if(2 * nd ~= rowB)
    error('MESS:error_arguments', 'number of rows of B differs from number of cols of A ( 2 * nd)');
end
%% solve
if isSym
    
    switch opB
        
        
        case 'N'
            X2 = eqn.K_ \ [B(nd + 1 : end, :); zeros(na, colB)];
            X = [B(1 : nd, :); - X2(1 : nd, :)];
            
            
        case 'T'
            X2 = eqn.K_ \ [B( : ,nd + 1 : end)'; zeros(na, colB)];
            X = [B( : , 1 : nd)'; - X2(1 : nd, :)];
    end
    
else
    switch opA
        
        case 'N'
            switch opB
                
                
                case 'N'
                    X2 = eqn.K_ \ [B(nd + 1 : end, :); zeros(na, colB)];
                    X = [B(1 : nd, :); - X2(1 : nd, :)];
                    
                    
                case 'T'
                    X2 = eqn.K_ \ [B( : ,nd + 1 : end)'; zeros(na, colB)];
                    X = [B( : , 1 : nd)'; - X2(1 : nd, :)];
            end
            
        case 'T'
            switch opB
                
                
                case 'N'
                    X2 = eqn.K_' \ [B(nd + 1 : end, :); zeros(na, colB)];
                    X = [B(1 : nd, :); - X2(1 : nd, :)];
                    
                    
                case 'T'
                    X2 = eqn.K_' \ [B( : ,nd + 1 : end)'; zeros(na, colB)];
                    X = [B( : , 1 : nd)'; - X2(1 : nd, :)];
            end
            
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Performace tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   for opA = 'N' compared three variants:
%   Var 1: ( Schur-Complement)
%                 X2 = (eqn.K_(1 : nd, 1 : nd) + eqn.K_(1 : nd, nd + 1 : end) * ...
%                     (eqn.K_(nd + 1 : end, nd + 1 : end) \ eqn.K_(nd + 1 : end, 1 : nd)))...
%                     \ B(nd + 1 : end, : );
%                 X = [B(1 : nd, :); - X2];
%   Var 2: (Schur-Complement + K_na_nd = K_nd_na^T)
%                 X2 = (eqn.K_(1 : nd, 1 : nd) + eqn.K_(1 : nd, nd + 1 : end) * ...
%                     (eqn.K_(nd + 1 : end, nd + 1 : end) \ eqn.K_(1 : nd, nd + 1 : end)'))...
%                     \ B(nd + 1 : end, : );
%                 X = [B(1 : nd, :); - X2];
%   Var 3: (direct exploit sparsity of K_)
%                 X2 = eqn.K_ \ [B(nd + 1 : end, :); zeros(na, cols)];
%                 X = [B(1 : nd, :); - X2(1 : nd, :)];
%   --> Variant 3 fastest (Jul. 2013 MATLAB R2012a 7.14.0.739)
