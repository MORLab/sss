function C=mul_ApE_default(eqn, opts,opA,p,opE,B,opB)

% function C=mul_ApE_default(eqn, opts,opA,p,opE,B,opB)
%
% This function returns C = (A_+p*E_)*B, where matrices A_ and E_ given by a structure eqn and input matrix B could be transposed.
%
%   Inputs:
%
%   eqn     structure containing fields 'A_' and 'E_'
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A_
%           opA = 'N' performs (A_ + p*opE(E_))*opB(B)
%           opA = 'T' performs (A_' + p*opE(E_))*opB(B)
%   p       scalar value
%   opE     character specifying the shape of E_
%           opA = 'N' performs (opA(A_) + p*E_)*opB(B)
%           opA = 'T' performs (opA(A_) + p*E_')*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs (opA(A_) + p*opE(E_))*B
%           opB = 'T' performs (opA(A_) + p*opE(E_))*B'
%
%   Output:
%
%   C = (opA(A_)+p*opE(E_))*opB(B)
%
% This function uses another default function mul_A_default(eqn,opA,B,opB) to obtain the result if E=I.

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
if (~ischar(opA) || ~ischar(opB) || ~ischar(opE))
    error('MESS:error_arguments', 'opA, opB or opE is not a char');
end

opA = upper(opA); opB = upper(opB); opE = upper(opE);

if(~(opA=='N' || opA=='T'))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(~(opB=='N' || opB=='T'))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if(~(opE=='N' || opE=='T'))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(~isnumeric(p)) || (length(p) ~= 1)
   error('MESS:error_arguments','p is not numeric'); 
end

if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end

if (~isnumeric(B)) || (~ismatrix(B))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(eqn.haveE ==1)
    if(~isfield(eqn,'E_') || ~isfield(eqn,'A_'))
        error('MESS:error_arguments','field eqn.E_ or eqn.A_ is not defined');
    end
else
    if(~isfield(eqn,'A_'))
        error('MESS:error_arguments','field eqn.A_ is not defined');
    end    
end

[rowA,colA]=size(eqn.A_);

if(eqn.haveE==1)
%% perform solve operations for E_ ~= Identity
    switch opA
        
        case 'N'
            switch opE
                
                case 'N'
                    
                    switch opB
                        
                        %implement multiplication (A_+p*E_)*B=C
                        case 'N'
                            
                            if (colA ~= size(B,1))
                                error('MESS:error_arguments','number of columns of A_ differs with number of rows of B'); 
                            end
                            
                            C = (eqn.A_+p*eqn.E_)*B;
                           
                        %implement multiplication (A_+p*E_)*B'=C
                        case 'T'
                            
                            if(colA~=size(B,2))
                                error('MESS:error_arguments','number of columns of A_ differs with number of columns of B');
                            end
                            
                            C = (eqn.A_+p*eqn.E_)*B';
                            
                    end
                    
                case 'T'
                    
                    switch opB
                        
                        %implement multiplication (A_+p*E_')*B=C
                        case 'N'
                            
                            if (colA ~= size(B,1))
                                error('MESS:error_arguments','number of columns of A_ differs with number of rows of B'); 
                            end
                            
                            C = (eqn.A_+p*eqn.E_')*B;
                           
                        %implement multiplication (A_+p*E_')*B'=C
                        case 'T'
                            
                            if(colA~=size(B,2))
                                error('MESS:error_arguments','number of columns of A_ differs with number of columns of B');
                            end
                            
                            C = (eqn.A_+p*eqn.E_')*B';
                            
                    end
                    
            end
            
        case 'T'
            switch opE
                
                case 'N'
                    
                    switch opB
                        
                        %implement multiplication (A_'+p*E_)*B=C
                        case 'N'
                            
                            if (rowA ~= size(B,1))
                                error('MESS:error_arguments','number of rows of A_ differs with number of rows of B'); 
                            end
                            
                            C = (eqn.A_'+p*eqn.E_)*B;
                           
                        %implement multiplication (A_'+p*E_)*B'=C
                        case 'T'
                            
                            if (rowA ~= size(B,2))
                                error('MESS:error_arguments','number of rows of A_ differs with number of columns of B'); 
                            end
                            
                            C = (eqn.A_'+p*eqn.E_)*B';
                            
                    end
                    
                case 'T'
                    
                    switch opB
                        
                        %implement multiplication (A_'+p*E_')*B=C
                        case 'N'
                            
                            if (rowA ~= size(B,1))
                                error('MESS:error_arguments','number of rows of A_ differs with number of rows of B'); 
                            end
                            
                            C = (eqn.A_'+p*eqn.E_')*B;
                           
                        %implement multiplication (A_'+p*E_')*B'=C
                        case 'T'
                            
                            if (rowA ~= size(B,2))
                                error('MESS:error_arguments','number of rows of A_ differs with number of columns of B'); 
                            end
                            
                            C = (eqn.A_'+p*eqn.E_')*B';
                            
                    end
                    
            end
    end
elseif(eqn.haveE==0)
%% perform solve operations for E_ = Identity
    switch opA
        
        case 'N'
            
            switch opB
                
                %implement solve (A_+p*E_)*B=C
                case 'N'
                    
                    if (colA ~= size(B,1))
                        error('MESS:error_arguments','number of columns of A_ differs with number of rows of B'); 
                    end
                    
                    C = mul_A_default(eqn,opts,'N',B,'N')+p*B;
                    
                %implement solve (A_+p*E_)*B'=C
                case 'T'
                    
                    if (colA ~= size(B,2))
                        error('MESS:error_arguments','number of columns of A_ differs with number of columns of B'); 
                    end
                    
                    C = mul_A_default(eqn,opts,'N',B,'T')+p*B';
                    
            end
            
        case 'T'
            
            switch opB
                
                %implement solve (A_'+p*E_)*B=C
                case 'N'
                    
                    if (rowA ~= size(B,1))
                        error('MESS:error_arguments','number of rows of A_ differs with number of rows of B'); 
                    end
                    
                    C = mul_A_default(eqn,opts,'T',B,'N')+p*B;
                    
                %implement solve (A_'+p*E_)*B'=C
                case 'T'
                    
                    if (rowA ~= size(B,2))
                        error('MESS:error_arguments','number of rows of A_ differs with number of columns of B'); 
                    end
                    
                    C = mul_A_default(eqn,opts,'T',B,'T')+p*B';
                    
            end
            
    end
end
end


