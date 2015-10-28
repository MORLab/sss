function [data,X,tx] = sim(sys,data,method,Ts_sample)
% Simulates a sss system using iddata input time series
% ------------------------------------------------------------------
% [data,x] = sim(sys,data)
% Inputs:       * sys: an sss-object containing the LTI system
%               * data: iddata object containing input time series
%               * method: time ingetration method (available are:
%               'forwardEuler', 'backwardEuler', 'RK4', and 'discrete'
%                default: continues 'RK4', discrete 'discrete'
%               * Ts_sample: sampling rate of the state-vector
% Outputs:      * data: iddata object containing output time series
%               * X: matrix of state vectors
%               * tx: time vector for X
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Stefan Jaensch (jaensch@tfd.mw.tum.de)
% Last Change:
% ------------------------------------------------------------------

if ~any(cellfun(@(x) strcmp('',x),sys.u))
    sys = sys.truncate(':',data.InputName);
else
    warning('InputName of sys is not specifing. Ignoring InputName definition in iddata')
end
[A,B,C,D,E] = ABCDE(sys);

x = sys.x0;
u = data.u;

Ts = data.Ts;
if isempty(Ts)
    error('nonconstant time steps are not supported');
end

if sys.Ts ~= 0 && sys.Ts~=Ts
    error('Sampling time of model and data must match')
end
if nargin==2 
    if sys.Ts ==0
        method = 'RK4';        
    else
        method = 'discrete';
    end      
    Ts_sample = Ts;
elseif nargin==3
    Ts_sample = Ts;
end
switch method
    case 'forwardEuler'
        if nargout == 1
            y = sss.sim_forwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        else
            [y,X,index] = sss.sim_forwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
            tx = data.SamplingInstants(index);
        end
    case 'backwardEuler'
        if nargout == 1
            y = sss.sim_backwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        else
            [y,X,index] = sss.sim_backwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
            tx = data.SamplingInstants(index);
        end
    case 'RK4'
        if nargout == 1
            y = sss.sim_RK4(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        else
            [y,X,index] = sss.sim_RK4(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
            tx = data.SamplingInstants(index);
        end
    case 'discrete'
        if nargout == 1
            y = sss.sim_discrete(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
        else
            [y,X,index] = sss.sim_discrete(A,B,C,D,E,u,x,Ts,Ts_sample,sys.isDescriptor);
            tx = data.SamplingInstants(index);
        end
end

data = [data iddata(y',[],Ts,'OutputName',sys.OutputName,'Tstart',data.Tstart)];
