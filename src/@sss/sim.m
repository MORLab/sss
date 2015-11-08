function [data,X,tx] = sim(sys,data,method,Ts_sample)
% sim - Simulates a sss system using iddata input time series
%
% Syntax:
%       [data,x,tx] = sim(sys,data)
%       [data,x,tx] = sim(sys,data,method)
%       [data,x,tx] = sim(sys,data,method,Ts_sample)
%
% Description:
%       Simulates a sss system using iddata input time series
%
% Input Arguments:       
%       *Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%       -data:  iddata object containing input time series
%       *Optional Input Arguments:*
%       -method:    time ingetration method (available are:
%                   'forwardEuler', 'backwardEuler', 'RK4', and 'discrete'
%                   default: continues 'RK4', discrete 'discrete'
%       -Ts_sample: sampling rate of the state-vector
%
% Output Arguments:      
%       -data: iddata object containing output time series
%       -X: matrix of state vectors
%       -tx: time vector for X
%
% Examples:
%       TODO
%
% See Also:
%       iddata/sim, DynamicSystem/lsim
%
% References:
%       TODO
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
% Authors:      Stefan Jaensch
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if ~any(cellfun(@(x) strcmp('',x),sys.u))
    sys = sys.truncate(':',data.InputName);
else
    warning('InputName of sys is not specifing. Ignoring InputName definition in iddata')
end
[A,B,C,D,E] = dssdata(sys);

x = sys.x0;
u = data.u;

Ts = data.Ts;
if isempty(Ts)
    error('nonconstant time steps are not supported');
end

if sys.Ts ~= 0 && sys.Ts ~= Ts
    error('Sampling time of model and data must match')
end
if nargin == 2 
    if sys.Ts == 0
        method = 'RK4';        
    else
        method = 'discrete';
    end      
    Ts_sample = Ts;
elseif nargin == 3
    Ts_sample = Ts;
end
switch method
    case 'forwardEuler'
        if nargout == 1
            y = sss.simForwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        else
            [y,X,index] = sss.simForwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
            tx = data.SamplingInstants(index);
        end
    case 'backwardEuler'
        if nargout == 1
            y = sss.simBackwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample);
        else
            [y,X,index] = sss.simBackwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample);
            tx = data.SamplingInstants(index);
        end
    case 'RK4'
        if nargout == 1
            y = sss.simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        else
            [y,X,index] = sss.simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
            tx = data.SamplingInstants(index);
        end
    case 'discrete'
        if nargout == 1
            y = sss.simDiscrete(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        else
            [y,X,index] = sss.simDiscrete(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
            tx = data.SamplingInstants(index);
        end
end

data = [data iddata(y',[],Ts,'OutputName',sys.OutputName,'Tstart',data.Tstart)];
