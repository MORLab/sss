function sys = c2d(sys,Ts,method)
% C2D - Converts a sss object from continues to discrete
%
% Syntax:
%       sys = C2D(sys,Ts,method)
%
% Description:
%       C2D coverts a continuous sparse state system sys in a discrete sparse
%       state system. It is required for the conversion the samping Time Ts
%       and the chosen method, that can be explicit Euler or implicit Euler.
%
% Input Arguments:
%       -sys: continuous time sss-object
%       -Ts:  sampling time
%       -method: string containing the selected discretization method.
%       Possible options are 'forward' (explicit Euler) or 'backward'
%       (implicit Euler)
%
% Output Arguments:
%       -sys: discrete time sss-object
%
% Examples:
%       The benchmark "build" is loaded and converted to a discrete model:
%
%> load build.mat
%> sysC = sss(A,B,C)
%> sysD = c2d(sysC,0.001,'forward')
%
% See Also:
%       ss/c2d
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Chair of Thermofluid Dynamics, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Stefan Jaensch (jaensch@tfd.mw.tum.de)
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  04 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if sys.Ts ~= 0
    error('Only continues models can be transformed')
end 
   
if nargin==2
    method = 'forward'; %default method: 'forward'
end
switch method
    case 'forward' % s = (z-1)/Ts
        sys.A = sys.E + Ts * sys.A;
        sys.B =  Ts * sys.B;        
    case 'backward' % s = (z-1)/(Ts*z)
        A = sys.A; E = sys.E;
        sys.E = E - Ts * A;
        sys.A = E;
        sys.B =  Ts * sys.B;
    otherwise
        error(['Method: ' method ' is not defined'])
end
sys.Ts = Ts;
