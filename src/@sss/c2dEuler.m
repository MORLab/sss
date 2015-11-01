function sys = c2dEuler(sys,Ts,method)
% Converts a sss object from continues to discrete using forward Euler
% scheme
% ------------------------------------------------------------------
% [sys] = c2d_Euler(sys, Ts,method)
% Inputs:       * sys: continuous time sss-object
%               * Ts: Sampling time
%               * method: forward = explicit Euler (default)
%                         backward = implicit Euler (default)
% Outputs:      * sys: discrete time sss-object
% ------------------------------------------------------------------
% Authors:      Stefan Jaensch
%
% see also: freqresp
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
