function sys = c2d(sys,Ts,method)
% C2D - Converts a sss object from continues to discrete
%
% Syntax:
%       sys = C2D(sys,Ts)
%       sys = C2D(sys,Ts,method)
%
% Description:
%       C2D coverts a continuous sparse state system sys in a discrete sparse
%       state system. For the conversion, the samping time |Ts| is
%       required.
%
%       If the function is called without giving a method, then sys = c2d(sys,Ts) 
%       uses explicit Euler ('forward') as a default discretization method.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: continuous time sss-object
%       -Ts:  sampling time
%       *Optional Input Arguments:*
%       -method: string containing the selected discretization method. 
%                Possible options are: [{'forward'} / 'backward' / 'tustin' / 'zoh']
%
%//Note: This function currently works only for the discretization method 'forward',
%       but the other discretization methods will be implemented in one of the next
%       releases.
%
% Output Arguments:
%       -sys: discrete time sss-object
%
% Examples:
%       The benchmark 'building' is loaded and converted to a discrete model:
%
%> load building.mat
%> sysC = sss(A,B,C)
%> sysD = c2d(sysC,0.001,'forward')
%
% See Also:
%       ss/c2d
%
% References:
%       * *[1] Franklin, Powell, and Workman (1997)*, Digital Control of Dynamic Systems (3rd Edition), Prentice Hall.
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Professur fuer Thermofluiddynamik, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Stefan Jaensch (jaensch@tfd.mw.tum.de), Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
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
        %TODO
        error('The conversion with the discretization method backward is not implemented yet');
%         A = sys.A; E = sys.E;
%         sys.E = E - Ts * A;
%         sys.A = E;
%         sys.B =  Ts * sys.B;
    case 'tustin'
        %TODO
        error('The conversion with the discretization method tustin is not implemented yet');
    case 'zoh'
        %TODO
        error('The conversion with the discretization method zoh is not implemented yet');
    otherwise
        error(['Method: ' method ' is not defined'])
end
sys.Ts = Ts;
