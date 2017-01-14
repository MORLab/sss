function X=sol_ApE_default(eqn, opts,opA,p,opE,C,opC)

% function X=sol_ApE_default(eqn, opts,opA,p,opE,C,opC)
%
% This function returns X = (A_ + p*E_)\C, where matrices A_ and E_ given by structure eqn and input matrix C could be transposed.
% Matrices A_ and E_ are assumed to be quadratic.
%
%   Inputs:
%
%   eqn     structure containing fields 'A_' and 'E_'
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' solves (A_ + p* opE(E_))*X = opC(C) 
%           opA = 'T' solves (A_' + p* opE(E_))*X = opC(C) 
%   p       scalar value
%   opE     character specifying the shape of E_
%           opE = 'N' solves (opA(A_) + p* E_)*X = opC(C) 
%           opE = 'T' solves (opA(A_) + p* E_')*X = opC(C) 
%   C       n-x-p matrix
%   opC     character specifies the form of opC(C)
%           opC = 'N' solves (opA(A_) + p* opE(E_))*X = C
%           opC = 'T' solves (opA(A_) + p* opE(E_))*X = C'
%
%   Output:
%
%   X       matrix fullfilling equation (opA(A_)+p*opE(E_))*X = opC(C)
%
% This function does not use other default functions.

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
if (~ischar(opA) || ~ischar(opE) || ~ischar(opC))
    error('MESS:error_arguments', 'opA, opE or opC is not a char');
end

opA = upper(opA); opE = upper(opE); opC = upper(opC);

if(~(opA=='N' || opA=='T'))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(~(opE=='N' || opE=='T'))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(~(opC=='N' || opC=='T'))
    error('MESS:error_arguments','opC is not ''N'' or ''T''');
end

if(~isnumeric(p)) || (length(p) ~= 1)
    error('MESS:error_arguments','p is not numeric');
end

if (~isnumeric(C)) || (~ismatrix(C))
    error('MESS:error_arguments','C has to ba a matrix');
end

%% check data in eqn structure
if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
if(eqn.haveE ==1)
    if(~isfield(eqn,'E_') || ~isfield(eqn,'A_'))
        error('MESS:error_arguments','field eqn.E_ or eqn.A_ is not defined');
    end
else
    if(~isfield(eqn,'A_'))
        error('MESS:error_arguments','field eqn.A_ is not defined');
    end    
end

[rowA,colA] = size(eqn.A_);


if(eqn.haveE==1)
%% perfom solve operations for E_ ~= Identity
switch opA
    
    case 'N'
        
        switch opE
        
            case 'N'
          
                switch opC
            
                    %implement solve (A_+p*E_)*X=C
                    case 'N'
                        
                        if(rowA~=size(C,1)) 
                            error('MESS:error_arguments','number of rows of A differs with number of rows of C');
                        end
                
                        X = (eqn.A_+p*eqn.E_)\C;
              
                    %implement solve (A_+p*E_)*X=C'
                    case 'T'

                        if(rowA~=size(C,2)) 
                            error('MESS:error_arguments','number of rows of A differs with number of columns of C');
                        end
                        
                        X = (eqn.A_+p*eqn.E_)\C';
              
                end
                
            case 'T'
          
                switch opC
            
                    %implement solve (A_+p*E_')*X=C
                    case 'N'
                        
                        if(rowA~=size(C,1)) 
                            error('MESS:error_arguments','number of rows of A differs with number of rows of C');
                        end
                        
                        X = (eqn.A_+p*eqn.E_')\C;
              
                    %implement solve (A_+p*E_')*X=C'
                    case 'T'
                        
                        if(rowA~=size(C,2)) 
                            error('MESS:error_arguments','number of rows of A differs with number of columns of C');
                        end
                        
                        X = (eqn.A_+p*eqn.E_')\C';
              
                end
                
        end
        
    case 'T'
        
        switch opE
        
            case 'N'
          
                switch opC
            
                    %implement solve (A_'+p*E_)*X=C
                    case 'N'
                        
                        if(colA~=size(C,1)) 
                            error('MESS:error_arguments','number of columns of A differs with number of rows of C');
                        end
                        
                        X = (eqn.A_'+p*eqn.E_)\C;
             
                    %implement solve (A_'+p*E_)*X=C'
                    case 'T'
                        
                        if(colA~=size(C,2)) 
                            error('MESS:error_arguments','number of columns of A differs with number of columns of C');
                        end
                        
                        X = (eqn.A_'+p*eqn.E_)\C';
              
                end
                
            case 'T'
          
                switch opC
            
                    %implement solve (A_'+p*E_')*X=C
                    case 'N'
                        
                        if(colA~=size(C,1)) 
                            error('MESS:error_arguments','number of columns of A differs with number of rows of C');
                        end
                        
                        X = (eqn.A_'+p*eqn.E_')\C;
              
                    %implement solve (A_'+p*E_')*X=C'
                    case 'T'
                        
                        if(colA~=size(C,2)) 
                            error('MESS:error_arguments','number of columns of A differs with number of columns of C');
                        end
                        
                        X = (eqn.A_'+p*eqn.E_')\C';
              
                end
        end
        
end
elseif(eqn.haveE==0)
  %% perform solve operations for E_ = Identity
  %% Note that init has put E_=speye 
  %% Note further that E_ is never transposed in contrast to the above
  
  switch opA
    
      case 'N'
      
          switch opC
        
              %implement solve (A_+p*I)*X=C
              case 'N'
                  
                  if(rowA~=size(C,1)) 
                      error('MESS:error_arguments','number of rows of A differs with number of rows of C');
                  end
                  
                  X = (eqn.A_+p*eqn.E_)\C;
          
              %implement solve (A_+p*I)*X=C'
              case 'T'
                  
                  if(rowA~=size(C,2)) 
                      error('MESS:error_arguments','number of rows of A differs with number of columns of C');
                  end
                  
                  X = (eqn.A_+p*eqn.E_)\C';
          
          end
          
      case 'T'
      
          switch opC
        
              %implement solve (A_'+p*I)*X=C
              case 'N'
                  
                  if(colA~=size(C,1)) 
                      error('MESS:error_arguments','number of columns of A differs with number of rows of C');
                  end
                  
                  X = (eqn.A_'+p*eqn.E_)\C;
          
              %implement solve (A_'+p*I)*X=C'
              case 'T'
                  
                  if(colA~=size(C,2)) 
                      error('MESS:error_arguments','number of columns of A differs with number of columns of C');
                  end
                  
                  X = (eqn.A_'+p*eqn.E_)\C';
          
          end
          
  end
end
end
