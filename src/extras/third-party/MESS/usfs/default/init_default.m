function [eqn, erg] = init_default(eqn, opts,flag1,flag2)

% function [eqn, erg] = init_default(eqn, opts,flag1,flag2)
%
% The function returns true or false if data for A_ and E_ resp. flag1 and flag2  are availabe and corrects in structure eqn.
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
%   This function does not use other default functions.
%
%   This function calls two other functions checkA and checkE implemented at the end.
%
%   The function checkA(eqn) proofs if a field 'A_' is included in the structure eqn and if the field 'A_' is numeric and quadratic.
%   This function returns the changed structure eqn and a boolean value erg (1- 'A_' is in structure eqn and a numeric and quadratic field)
%
%   The function checkE(eqn) proofs if a field 'E_' is included in the structure eqn and if the field 'E_' is numeric and quadratic.
%   If the structure does not include a field E, a new field 'E_' is defined as a sparse identity matrix by size of field 'A_'.
%   This function returns the changed structure eqn and a boolean value erg (1- 'E_' is in structure eqn and a numeric and quadratic field)

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

%erg = init_default(eqn, flag1);    
elseif(na==3)
    switch flag1
        case {'A','a'}
            [eqn,erg] = checkA(eqn);
        case {'E','e'}
            [eqn,erg] = checkE(eqn);
        otherwise
            error('MESS:control_data','flag1 has to be ''A_'' or ''E_''');
    end
    
%erg = init_default(eqn,flag1,flag2);
elseif(na==4)
    switch flag1
        case {'A','a'}
            [eqn, erg] = checkA(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA] = checkA(eqn);
                    erg = erg && ergA;
                case {'E','e'}
                    [eqn, ergE]= checkE(eqn);
                    erg = erg && ergE; 
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        case {'E','e'}
            [eqn, erg] = checkE(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA] = checkA(eqn);
                    erg = erg && ergA;
                case {'E','e'}
                    [eqn, ergE]= checkE(eqn);
                    erg = erg && ergE; 
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end 
end

end

%checkdata for A_
function [eqn, erg] = checkA(eqn)
erg = isfield(eqn,'A_');
if(erg)
    erg = isnumeric(eqn.A_);
end
erg=erg&&(size(eqn.A_,1)==size(eqn.A_,2));
end

%checkdata for E_
function [eqn, erg] = checkE(eqn)
if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
if ~eqn.haveE
  erg = 1;
  eqn.E_= speye(size(eqn.A_,1)); %make sure we have an identity for
                                 %computations in ApE functions
else
  erg = isfield(eqn,'E_');
  if(erg)
    erg = isnumeric(eqn.E_);
  end
  erg=erg&&(size(eqn.E_,1)==size(eqn.E_,2));
end
end
