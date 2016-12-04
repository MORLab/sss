function [eqn,erg] = init_so_2(eqn, opts,flag1,flag2)

%function [eqn,erg] = init_so_2(eqn, opts,flag1,flag2)
%
% The second order system
%
%   M x"(t) + D x'(t) + K x(t) = B u(t)
%                         y(t) = C x(t)
%       
% is transformed to the first order system
%
%   E z'(t) = A z(t) + G u(t)
%      y(t) = L z(t)
%  
% where
%
%      | D  M|
%   E= | M  0|
%   
%      |-K  0|
%   A= | 0  M|
%   
%      | B |
%   G= | 0 |
%   
%   L= [C  0]
%   
%         | x(t)  |
%   z(t)= | x'(t) | .
%   
% Matrices M, D, K are assumed to be quadratic, symmetric and positive definit.
% The function returns true or false if data for A and E resp. flag1 and flag2  are availabe and corrects in structure eqn.
%
%   Inputs:
%
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   flag1           'A'/'E' to check if A or E is in eqn
%   flag2           'A'/'E' to check if A or E is in eqn
%
%   Outputs:
%
%   eqn             changed structure with data
%   erg             1 if data corresponding to flag1 (and flag2) are available , 0 data are not available 
%
%   This function does not use other so3 functions.
%
%   The function checkA(eqn) proofs if the fields 'K_' and 'M_' are included in the structure eqn and if these fields are numeric.
%   This function returns the changed structure eqn and a boolean value erg (1- 'K_' and 'M_' are in structure eqn and a numeric field)
%
%   The function checkE(eqn) proofs if the fields 'D_' and 'M_' are included in the structure eqn and if these fields are numeric.
%   This function returns the changed structure eqn and a boolean value erg (1-  'D_' and 'M_' are in structure eqn and a numeric field)

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

%start checking
na = nargin;
if(na<=2)
    error('MESS:control_data','Number of input Arguments are at least 3');

%erg = init_so_1(eqn, flag1);    
elseif(na==3)
    switch flag1
        case {'A','a'}
            [eqn,erg] = checkA(eqn);
        case {'E','e'}
            [eqn,erg] = checkE(eqn);
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end
    
%erg = init_so_1(eqn,flag1,flag2);
elseif(na==4)
    switch flag1
        case {'A','a'}
            [eqn,erg] = checkA(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA] = checkA(eqn);
                    erg = erg && ergA;
                case {'E','e'}
                    [eqn,ergE] = checkE(eqn);
                    erg = erg && ergE;
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        case {'E','e'}
            [eqn,erg] = checkE(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA] = checkA(eqn);
                    erg = erg &&ergA;
                case {'E','e'}
                    [eqn,ergE] = checkE(eqn);
                    erg =  erg && ergE;
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end 
end
end

%checkdata for A
function [eqn,erg] = checkA(eqn)
erg = isfield(eqn,'K_') &&isfield(eqn,'M_');
if(erg)
    erg = isnumeric(eqn.K_) && isnumeric(eqn.M_);
    if(~issparse(eqn.K_))
        warning('MESS:control_data','K is not sparse');
    end
    if(~issparse(eqn.M_))
        warning('MESS:control_data','M is not sparse');
    end
end
end

%checkdata for E
function [eqn,erg] = checkE(eqn)
if ~isfield(eqn, 'haveE'), eqn.haveE = 1; end
if ~eqn.haveE
    warning('MESS:control_data','eqn.haveE has to be 1');
end
erg = isfield(eqn,'D_')&&isfield(eqn,'M_');
if(erg)
    erg = isnumeric(eqn.D_) && isnumeric(eqn.M_);
    if(~issparse(eqn.M_))
        warning('MESS:control_data','M is not sparse');
    end
    if(~issparse(eqn.D_))
        warning('MESS:control_data','D is not sparse');
    end
end
end
