function [y,X,tx] = lsim(sys,u,Ts,method,Ts_sample)
% lsim - Simulate time response of sparse dynamic system to arbitrary inputs
%
% Syntax:
%       [y,X,tx] = lsim(sys,u,Ts)
%       [y,X,tx] = lsim(sys,u,Ts,Ts_sample)
%       [y,X,tx] = lsim(sys,u,Ts,Ts_sample,method)
%
% Description:
%       Simulates the response in time of sparse LTI systems to arbitrary inputs.
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%       -u: Input signal as vector
%       -Ts: Sampling rate of input vector
%       *Optional Input Arguments:*
%       -method:    time ingetration method
%                   [{'RK4' (continuous) / 'discrete' (discrete)} /
%                   'forwardEuler' / 'backwardEuler' / 'RKDP']
%       -Ts_sample: sampling rate of the state-vector
%
% Output Arguments:
%       -y: Matrix with the response of the system
%       -X: Matrix of state vectors
%       -tx: time vector for X
%
% See Also:
%       sss/sim
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Professur fuer Thermofluiddynamik, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Michael Leipold, Stefan Jaensch
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  17 Dec 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


[A,B,C,D,E] = dssdata(sys);
x = sys.x0;

if nargin == 2
    if sys.Ts == 0
        method = 'RK4';
    else
        method = 'discrete';
    end
    Ts_sample = Ts;
elseif nargin == 3
    Ts_sample = Ts;
elseif nargin == 4
    Ts_sample = inf;
end

switch method
    case 'forwardEuler'
        if nargout == 1
            y = simForwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        else
            [y,X,tx] = simForwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        end
    case 'backwardEuler'
        if nargout == 1
            y = simBackwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);            
        else
            [y,X,tx] = simBackwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        end
    case 'RK4'
        if nargout == 1
            y = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        else
            [y,X,tx] = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        end
    case 'discrete'
        if nargout == 1
            y = simDiscrete(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        else
            [y,X,tx] = simDiscrete(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        end
end

end

