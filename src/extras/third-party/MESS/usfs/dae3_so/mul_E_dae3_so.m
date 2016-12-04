function C = mul_E_dae3_so(eqn, opts, opE, B, opB)
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
for mat='MG'
    if(~isfield(eqn, sprintf('%c_',mat)) || ~eval(sprintf('isnumeric(eqn.%c_)',mat)))
        error('MESS:error_arguments', 'field eqn.%c_ is not defined',mat);
    end
end

%% perfom multiplication
nv = size(eqn.M_,1);
np = size(eqn.G_,1);

if (opB=='N' && (size(B,1)==(2*nv+np))) || (opB=='T' && (size(B,2)==(2*nv+np)))
    switch opE
        
        case 'N'
            switch opB
                
                %implement operation E_*B
                case 'N'
                    C = [B(1:nv,:);
                        eqn.M_*B(nv+1:2*nv,:)+eqn.alpha*eqn.G_'*B(2*nv+1:end,:);
                        eqn.alpha*eqn.G_*B(nv+1:2*nv,:)];
                    
                    %implement operation E_*B'
                case 'T'
                    C = [B(:,1:nv)';
                        eqn.M_*B(:,nv+1:2*nv)'+eqn.alpha*eqn.G_'*B(:,2*nv+1:end)';
                        eqn.alpha*eqn.G_*B(:,nv+1:2*nv)'];
                    
            end
            
        case 'T'
            
            switch opB
                %implement operation E_'*B
                case 'N'
                    C = [B(1:nv,:);
                        eqn.M_'*B(nv+1:2*nv,:)+eqn.alpha*eqn.G_'*B(2*nv+1:end,:);
                        eqn.alpha*eqn.G_*B(nv+1:2*nv,:)];
                    
                    %implement operation E_'*B'
                case 'T'
                    C = [B(:,1:nv)';
                        eqn.M_'*B(:,nv+1:2*nv)'+eqn.alpha*eqn.G_'*B(:,2*nv+1:end)';
                        eqn.alpha*eqn.G_*B(:,nv+1:2*nv)'];
                    
            end
            
    end
    
elseif (opB=='N' && (size(B,1)==(2*nv))) || (opB=='T' && (size(B,2)==(2*nv)))
    switch opE
        
        case 'N'
            switch opB
                
                %implement operation E_*B
                case 'N'
                    C = [B(1:nv,:);
                        eqn.M_*B(nv+1:2*nv,:)];
                    
                    %implement operation E_*B'
                case 'T'
                    C = [B(:,1:nv)';
                        eqn.M_*B(:,nv+1:2*nv)'];
                    
            end
            
        case 'T'
            
            switch opB
                %implement operation E_'*B
                case 'N'
                    C = [B(1:nv,:);
                        eqn.M_'*B(nv+1:2*nv,:)];
                    
                    %implement operation E_'*B'
                case 'T'
                    C = [B(:,1:nv)';
                        eqn.M_'*B(:,nv+1:2*nv)'];
            end
            
    end
else
    error('MESS:error_arguemnts', 'B has wrong number of cols');
end

end
