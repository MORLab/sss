function C=mul_E_dae_so_1(eqn, opts, opE, B, opB)

%% function mul_A_so_1 perfoms operation C = opE(E)*opB(B)
%
% Input:
%   eqn     structure contains field E
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' performs E*opB(B)
%           opE = 'T' performs E'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opE(E)*B
%           opB = 'T' performs opE(E)*B'
%
% Output:
% C = opE(E)*opB(B)
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
if (~isfield(eqn,'M_') || ~isnumeric(eqn.M_))
    error('MESS:equation_data',...
        'Empty or Corrupted field M detected in equation structure.')
elseif (~isfield(eqn,'D_') || ~isnumeric(eqn.D_))
    error('MESS:equation_data',...
        'Empty or Corrupted field D detected in equation structure.')
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
    [rowB, colB] = size(B);
else
    [colB, rowB] = size(B);
end

if(2 * nd ~= rowB)
    error('MESS:error_arguments', 'number of rows of B differs from number of cols of E ( 2 * nd)');
end


%% perfom multiplication
if isSym
    switch opE
        
        case 'N'
            switch opB
                
                
                case 'N'

                    C1 = B(nd + 1 : end, :);
                    C2 = eqn.M_(1:nd,1:nd) * B(1:nd,:) ...
                       + eqn.D_(1:nd,1:nd) * B(nd + 1 : end, :);
                    C = [C1; C2];
                    
                    
                case 'T'
                    C1 = B( : , nd + 1 : end)';
                    C2 = eqn.M_(1:nd,1:nd) * B(:,1:nd)' ...
                       + eqn.D_(1:nd,1:nd) * B( : , nd + 1 : end)';
                    C = [C1; C2];
            end
            
        case 'T'
            switch opB
                
                
                case 'N'
                    C1 = eqn.M_(1:nd,1:nd) * B(nd + 1 : end, :);
                    C2 = eqn.D_(1:nd,1:nd) * B(nd + 1 : end, :);
                    C = [C1; B(1 : nd, :) + C2];
                    
                    
                case 'T'
                    C1 = eqn.M_(1:nd,1:nd) * B( : , nd + 1 : end)';
                    C2 = eqn.D_(1:nd,1:nd) * B( : , nd + 1 : end)';
                    C = [C1; B( : , 1 : nd)' + C2];
            end
            
    end
else
    switch opE
        
        case 'N'
            switch opB
                                
                case 'N'

                    C1 = B(nd + 1 : end, :);
                    C2 = eqn.M_(1:nd,1:nd) * B(1:nd ,:) ...
                       + eqn.D_(1:nd,1:nd) * B(nd + 1 : end, :);
                    C = [C1; C2];
                    
                    
                case 'T'
                    C1 = B( : , nd + 1 : end)';
                    C2 = eqn.M_(1:nd,1:nd) * B(:,1:nd)'...
                       + eqn.D_(1:nd,1:nd) * B( : , nd + 1 : end)';
                    C = [C1; C2(1 : nd, :)];
            end
            
        case 'T'
            switch opB
                
                case 'N'
                    C1 = eqn.M_(1:nd,1:nd)' * B(nd + 1 : end, :);
                    C2 = eqn.D_(1:nd,1:nd)' * B(nd + 1 : end, :);
                    C = [C1; B(1 : nd, :) + C2];
                    
                    
                case 'T'
                    C1 = eqn.M_(1:nd,1:nd)' * B( : , nd + 1 : end)';
                    C2 = eqn.D_(1:nd,1:nd)' * B( : , nd + 1 : end)';
                    C = [C1; B( : , 1 : nd)' + C2];
            end
            
    end
end

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Performace tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   for opE = 'N' tested four variants:
%   Var 1:
%                 C1 = V(nd+1 : end, :);
%                 C2 = eqn.M_(1 : nd, 1 : nd) * V(1 : nd, :) + eqn.D_(1 : nd, 1 : nd) * V(nd + 1 : end, :);
%                 C = [C1; C2];
%   Var 2:
%                 C = [sparse(nd,nd) speye(nd,nd);
%                     eqn.M_(1 : nd, 1 : nd) eqn.D_(1 : nd, 1 : nd)] * V;
%   Var 3:
%                 C1 = V(nd+1 : end, :);
%                 C2 = [eqn.M_(1 : nd, 1 : nd) eqn.D_(1 : nd, 1 : nd)] * V;
%                 C = [C1; C2(1 : nd, :)];
%   Var 4:
%                 C1 = V(nd+1 : end, :);
%                 C2 = eqn.M_ * [V; zeros(na - nd, cols)] + eqn.D_ * [V(nd + 1 : end, :); zeros(na, cols)];
%                 C = [C1; C2(1 : nd, :)];
%   --> Variant 4 fastest (Jul. 2013 MATLAB R2012a 7.14.0.739)
%
%   for opE = 'T' tested three variants:
%   Var 1:
%                 C1 = eqn.M_(1 : nd, 1 : nd) * V(nd+1 : end, :);
%                 C2 =  V(1 : nd, :) + eqn.D_(1 : nd, 1 : nd) * V(nd + 1 : end, :);
%                 C = [C1; C2];
%   Var 2:
%                 V2 = [V(nd + 1 : end, :); zeros(na, cols)];
%                 C1 = eqn.M_ * V2;
%                 C2 =  [V(1 : nd, :); zeros(na, cols)] + eqn.D_ * V2;
%                 C = [C1(1 : nd, :); C2(1 : nd, :)];
%   Var 3:
%                 V2 = [V(nd + 1 : end, :); zeros(na, cols)];
%                 C1 = eqn.M_ * V2;
%                 C2 =  eqn.D_ * V2;
%                 C = [C1(1 : nd, :); V(1 : nd, :) + C2(1 : nd, :)];
%   --> Variant 3 fastest (Jul. 2013 MATLAB R2012a 7.14.0.739)
