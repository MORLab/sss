function X = sol_E_dae2_so(eqn, opts, opE, B, opB)
%% function sol_E_dae2_so solves opE(E)*X = opB(B) resp. performs X=opE(E)\opB(B)
%
% Input:
%   eqn     structure contains data for E (M_)
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' solves E *X = opB(B)
%           opE = 'T' sovles E'*X = opB(B)
%
%   B       p-x-q matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' solves opE(E)*X = B
%           opB = 'T' solves opE(E)*X = B'
%
% Output
%
%   X       matrix fullfills equation opE(E)*X = opB(B)
%
%% check input Paramters

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
if(~isfield(eqn, 'M_')) || ~isnumeric(eqn.M_)
    error('MESS:error_arguments', 'field eqn.M_ is not defined');
end

nv = size(eqn.M_,1);
np = size(eqn.G_,1);

%% solve
if (opB=='N' && (size(B,1)==(2*nv+np))) || (opB=='T' && (size(B,2)==(2*nv+np)))
    
    switch opE
        
        case 'N'
            switch opB
                
                %implement solve E*X=B
                case 'N'
                    X = [B(1:nv,:);
                        [eqn.M_,eqn.alpha*eqn.G_'; ...
                        eqn.alpha*eqn.G_,sparse(np,np)] \ B(nv+1:end,:)
                        ];
                    
                    %implement solve A*X=B'
                case 'T'
                    X = [B(:,1:nv)';
                        [eqn.M_,eqn.alpha*eqn.G_'; ...
                        eqn.alpha*eqn.G_,sparse(np,np)] \ B(:,nv+1:end)'
                        ];
                    
            end
            
        case 'T'
            switch opB
                
                %implement solve E'*X=B
                case 'N'
                    X = [B(1:nv,:);
                        [eqn.M_',eqn.alpha*eqn.G_'; ...
                        eqn.alpha*eqn.G_,sparse(np,np)] \ B(nv+1:end,:)
                        ];
                    
                    %implement solve A_'*X=B'
                case 'T'
                    X = [B(:,1:nv)';
                        [eqn.M_',eqn.alpha*eqn.G_'; ...
                        eqn.alpha*eqn.G_,sparse(np,np)] \ B(:,nv+1:end)'
                        ];
            end
            
    end
    
elseif (opB=='N' && (size(B,1)==(2*nv))) || (opB=='T' && (size(B,2)==(2*nv)))
    error('MESS:error_usage','sol_E_dae2_so is only coded for shift parameter computation');
else
    error('MESS:error_arguemnts', 'B has wrong number of cols');
end


end
