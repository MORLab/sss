function C = mul_ApE_dae_1(eqn, opts, opA, p, opE, B, opB)

%% function mul_A perfoms operation C = (opA(A_)+pc*opE(E_))*opB(B)
%
% Input:
%   eqn     structure contains A_
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
% Output:
% C = (opA(A_)+pc*opE(E_))*opB(B)
%
%   uses size_dae_1

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
if (~ischar(opA) || ~ischar(opB) || ~ischar(opE))
    error('MESS:error_arguments', 'opA, opB or opE is not a char');
end

opA = upper(opA); opB = upper(opB); opE = upper(opE);

if(~(opA == 'N' || opA == 'T'))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(~(opB == 'N' || opB == 'T'))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if(~(opE == 'N' || opE == 'T'))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(~isnumeric(p))
   error('MESS:error_arguments','p is not numeric'); 
end

if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end

if (~isnumeric(B)) || (~ismatrix(B))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(eqn.haveE ==1)
    if(~isfield(eqn,'E_') || ~isnumeric(eqn.E_)...
            || ~isfield(eqn,'A_')) || ~isnumeric(eqn.A_)
        error('MESS:error_arguments','field eqn.E_ or eqn.A_ is not defined or corrupted');
    end
else
    if(~isfield(eqn,'A_')) || ~isnumeric(eqn.A_)
        error('MESS:error_arguments','field eqn.A_ is not defined');
    end    
end
if ~isfield(eqn, 'st')    || ~isnumeric(eqn.st)
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.')
end
n = size(eqn.A_, 1);
st = eqn.st;
[rowB,colB] = size(B);

if(opB == 'N')
    if(n > rowB)
        B = [B;
            zeros(n - st, colB)];
    elseif n < rowB
        error('MESS:error_arguments', 'B has more rows than A');
    end
else
    if(n > colB)
        B = [B, zeros(rowB, n - st)];
    elseif n < colB
        error('MESS:error_arguments', 'B has more columns than A');
    end
end

%% perfom multiplication for E_ ~= Identity
if(eqn.haveE == 1)
    switch opA
        
        case 'N'
            switch opE
                
                case 'N'
                    
                    switch opB
                        
                        %implement multiplication (A_+pc*E_)*B=C
                        case 'N'
                            C = (eqn.A_ + p * eqn.E_) * B;
                           
                        %implement multiplication (A_+pc*E_)*B'=C
                        case 'T'
                            C = (eqn.A_ + p * eqn.E_) * B';
                            
                    end
                    
                case 'T'
                    
                    switch opB
                        
                        %implement multiplication (A_+pc*E_')*B=C
                        case 'N'
                            C = (eqn.A_ + p * eqn.E_') * B;
                           
                        %implement multiplication (A_+pc*E_')*B'=C
                        case 'T'
                            C = (eqn.A_ + p * eqn.E_') * B';
                            
                    end
                    
            end
            
        case 'T'
            switch opE
                
                case 'N'
                    
                    switch opB
                        
                        %implement multiplication (A_+pc*E_)*B=C
                        case 'N'
                            C = (eqn.A_' + p * eqn.E_) * B;
                           
                        %implement multiplication (A_+pc*E_)*B'=C
                        case 'T'
                            C = (eqn.A_' + p * eqn.E_) * B';
                            
                    end
                    
                case 'T'
                    
                    switch opB
                        
                        %implement multiplication (A_+pc*E_')*B=C
                        case 'N'
                            C = (eqn.A_' + p * eqn.E_') * B;
                           
                        %implement multiplication (A_+pc*E_')*B'=C
                        case 'T'
                            C = (eqn.A_' + p * eqn.E_') * B';
                            
                    end
                    
            end
    end
elseif(eqn.haveE==0)
%% perform multiplication for E_ = Identity
    switch opA
        
        case 'N'
            
            switch opB
                
                %implement solve (A_+pc*E_)*B=C
                case 'N'
                    C = mul_A(eqn, 'N', B, 'N') + p * B;
                    
                %implement solve (A_+pc*E_)*B=C
                case 'T'
                    C = mul_A(eqn, 'N', B, 'T') + p * B';
                    
            end
            
        case 'T'
            
            switch opB
                
                %implement solve (A_'+pc*E_)*B=C
                case 'N'
                    C = mul_A(eqn, 'T', B, 'N') + p * B;
                    
                %implement solve (A_'+pc*E_)*B'=C
                case 'T'
                    C = mul_A(eqn, 'T', B, 'T') + p * B';
                    
            end
            
    end
end
C = C(1 : st, :);
end



