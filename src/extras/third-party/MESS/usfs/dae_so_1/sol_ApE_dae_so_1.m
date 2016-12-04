function X = sol_ApE_dae_so_1(eqn, opts, opA, p, opE, C, opC)

%% function sol_ApE_so_1 solves (opA(A) + p*opE(E))*X = opC(C) resp. performs X=(opA(A)+p*opE(E))\opC(C)
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
%   p      scalar Value
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' for E
%           opE = 'T' for E'
%
%   C       n-x-p matrix
%
%   opC     character specifies the form of opC(C)
%           opC = 'N' for C
%           opC = 'T' for C'
%
% Output
%
%   X       matrix fullfills equation (opA(A)+p*opE(E))*X = C
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
if (~ischar(opA) || ~ischar(opE) || ~ischar(opC))
    error('MESS:error_arguments', 'opA, opE or opC is not a char');
end

opA = upper(opA); opE = upper(opE); opC = upper(opC);

if(~(opA == 'N' || opA == 'T'))
    error('MESS:error_arguments', 'opA is not ''N'' or ''T''');
end

if(~(opE == 'N' || opE == 'T'))
    error('MESS:error_arguments', 'opE is not ''N'' or ''T''');
end

if(~(opC == 'N' || opC == 'T'))
    error('MESS:error_arguments', 'opC is not ''N'' or ''T''');
end
if(~isnumeric(p))
   error('MESS:error_arguments','p is not numeric'); 
end
if (~isnumeric(C)) || (~ismatrix(C))
    error('MESS:error_arguments','C has to ba a matrix');
end


%% check data in eqn structure
if (~isfield(eqn,'K_') || ~isnumeric(eqn.K_))
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.')
end
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
if ~isfield(eqn,'haveE'), eqn.haveE=0; end

[rowK, colK] = size(eqn.K_);
nd = eqn.nd;

if(opC == 'N')
    rowC = size(C, 1);
    colC = size(C, 2);
else
    rowC = size(C, 2);
    colC = size(C, 1);
end

if(2 * nd ~= rowC)
    error('MESS:error_arguments','Rows of A differs from rows of C');
end


%% solve (A + p * E) * x = C 
%% perfom solve operations for E ~= Identity
if isSym
    if(eqn.haveE == 1)
        switch opE
            
            case 'N'
                
                switch opC
                    
                    
                    case 'N'
                        X2 = (eqn.K_ + p * ( p * eqn.M_ - eqn.D_)) \ ...
                            (p * eqn.M_ * [C(1 : nd, :); zeros(rowK - nd, colC)]...
                            - [C(nd + 1 : end, :); zeros(rowK - nd, colC)]);
                        X2 = X2(1 : nd, :);
                        X1 = C(1 : nd, :) - p * X2;
                        X = [X1; X2];
                        
                        
                    case 'T'
                        X2 = (eqn.K_ + p * ( p * eqn.M_ - eqn.D_)) \ ...
                            (p * eqn.M_ * [C( : , 1 : nd)'; zeros(rowK - nd, colC)]...
                            - [C( : , nd + 1 : end)'; zeros(rowK - nd, colC)]);
                        X2 = X2(1 : nd, :);
                        X1 = C( : , 1 : nd)' - p * X2;
                        X = [X1; X2];
                        
                end
                
            case 'T'
                
                switch opC
                    
                    
                    case 'N'
                        X2 = (eqn.K_ + p * ( p * eqn.M_ - eqn.D_)) \ ...
                            (p * [C(1 : nd, :); zeros(rowK - nd, colC)]...
                            - [C(nd + 1 : end, :); zeros(rowK - nd, colC)]);
                        X2 = X2(1 : nd, :);
                        X1 = C(1 : nd, :) - p * eqn.M_(1 : nd, 1 : nd) * X2;
                        X = [X1; X2];
                        
                        
                    case 'T'
                        X2 = (eqn.K_ + p * ( p * eqn.M_ - eqn.D_)) \ ...
                            (p * [C( : , 1 : nd)'; zeros(rowK - nd, colC)]...
                            - [C( : , nd + 1 : end)'; zeros(rowK - nd, colC)]);
                        X2 = X2(1 : nd, :);
                        X1 = C( : , 1 : nd)' - p * eqn.M_(1 : nd, 1 : nd) * X2;
                        X = [X1; X2];
                        
                end
                    
        end
            
    elseif(eqn.haveE==0)
        %% perform solve operations for E_ = Identity
        
        
        switch opC
            
            
            case 'N'
                X1 = C(1 : nd, :) / (p + 1);
                X2 = (eqn.K_ - p * [speye(nd, nd) sparse(nd, colK - nd); ...
                    sparse(rowK - nd, colK)]) \ [-C(nd + 1 : end, :); zeros(rowK - nd, colC)];
                X2 = X2(1 : nd, :);
                X = [X1; X2];
                
                    
            case 'T'
                X1 = C( : , 1 : nd)' / (p + 1);
                X2 = (eqn.K_ - p * [speye(nd, nd) sparse(nd, colK - nd); ...
                    sparse(rowK - nd, colK)]) \ [-C( : , nd + 1 : end)'; zeros(rowK - nd, colC)];
                X2 = X2(1 : nd, :);
                X = [X1; X2];
                
        end
    end
else
    if(eqn.haveE == 1)
        switch opA
            
            case 'N'
                switch opE
                    
                    case 'N'
                        
                        switch opC
                            
                            
                            case 'N'
                                X2 = (eqn.K_ + p * ( p * eqn.M_ - eqn.D_)) \ ...
                                    (p * eqn.M_ * [C(1 : nd, :); zeros(rowK - nd, colC)]...
                                    - [C(nd + 1 : end, :); zeros(rowK - nd, colC)]);
                                X2 = X2(1 : nd, :);
                                X1 = C(1 : nd, :) - p * X2;
                                X = [X1; X2];
                                
                                
                            case 'T'
                                X2 = (eqn.K_ + p * ( p * eqn.M_ - eqn.D_)) \ ...
                                    (p * eqn.M_ * [C( : , 1 : nd)'; zeros(rowK - nd, colC)]...
                                    - [C( : , nd + 1 : end)'; zeros(rowK - nd, colC)]);
                                X2 = X2(1 : nd, :);
                                X1 = C( : , 1 : nd)' - p * X2;
                                X = [X1; X2];
                                
                        end
                            
                    case 'T'
                            
                        switch opC
                            
                            
                            case 'N'
                                X2 = (eqn.K_ + p * ( p * eqn.M_ - eqn.D_)') \ ...
                                    (p * [C(1 : nd, :); zeros(rowK - nd, colC)]...
                                    - [C(nd + 1 : end, :); zeros(rowK - nd, colC)]);
                                X2 = X2(1 : nd, :);
                                X1 = C(1 : nd, :) - p * eqn.M_(1 : nd, 1 : nd)' * X2;
                                X = [X1; X2];
                                
                                
                            case 'T'
                                X2 = (eqn.K_ + p * ( p * eqn.M_ - eqn.D_)') \ ...
                                    (p * [C( : , 1 : nd)'; zeros(rowK - nd, colC)]...
                                    - [C( : , nd + 1 : end)'; zeros(rowK - nd, colC)]);
                                X2 = X2(1 : nd, :);
                                X1 = C( : , 1 : nd)' - p * eqn.M_(1 : nd, 1 : nd)' * X2;
                                X = [X1; X2];
                                
                        end
                        
                end
                    
            case 'T'
                switch opE
                    
                    case 'N'
                            
                        switch opC
                            
                            
                            case 'N'
                                X2 = (eqn.K_' + p * ( p * eqn.M_ - eqn.D_)) \ ...
                                    (p * eqn.M_ * [C(1 : nd, :); zeros(rowK - nd, colC)]...
                                    - [C(nd + 1 : end, :); zeros(rowK - nd, colC)]);
                                X2 = X2(1 : nd, :);
                                X1 = C(1 : nd, :) - p * X2;
                                X = [X1; X2];
                                
                                
                            case 'T'
                                X2 = (eqn.K_' + p * ( p * eqn.M_ - eqn.D_)) \ ...
                                    (p * eqn.M_ * [C( : , 1 : nd)'; zeros(rowK - nd, colC)]...
                                    - [C( : , nd + 1 : end)'; zeros(rowK - nd, colC)]);
                                X2 = X2(1 : nd, :);
                                X1 = C( : , 1 : nd)' - p * X2;
                                X = [X1; X2];
                                
                        end
                            
                    case 'T'
                            
                        switch opC
                            
                                
                            case 'N'
                                X2 = (eqn.K_ + p * ( p * eqn.M_ - eqn.D_))' \ ...
                                    (p * [C(1 : nd, :); zeros(rowK - nd, colC)]...
                                    - [C(nd + 1 : end, :); zeros(rowK - nd, colC)]);
                                X2 = X2(1 : nd, :);
                                X1 = C(1 : nd, :) - p * eqn.M_(1 : nd, 1 : nd)' * X2;
                                X = [X1; X2];
                                
                                
                            case 'T'
                                X2 = (eqn.K_ + p * ( p * eqn.M_ - eqn.D_))' \ ...
                                    (p * [C( : , 1 : nd)'; zeros(rowK - nd, colC)]...
                                    - [C( : , nd + 1 : end)'; zeros(rowK - nd, colC)]);
                                X2 = X2(1 : nd, :);
                                X1 = C( : , 1 : nd)' - p *eqn.M_(1 : nd, 1 : nd)' * X2;
                                X = [X1; X2];
                                
                        end
                            
                end
                    
        end
   
    end
end
end


