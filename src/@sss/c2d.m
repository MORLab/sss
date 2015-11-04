function sys = c2d(sys,Ts,method)
% C2D - Converts a sss object from continues to discrete
%
% Syntax:
%       sys = C2D(sys,Ts,method)
%
% Description:
%       TODO
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
%       TODO
%
% See Also:
%       freqresp
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
    method = 'forward';
end
switch method
    case 'forward'
        sys.A = sys.E + Ts * sys.A;
        sys.B =  Ts * sys.B;        
    case 'backward'
        A = sys.A; E = sys.E;
        sys.E = E - Ts * A;
        sys.A = E;
        sys.B =  Ts * sys.B;
    otherwise
        error(['Method: ' method ' is not defined'])
end
sys.Ts = Ts;
