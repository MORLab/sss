function [data,X,tx] = sim(sys,data,method,Ts_sample)
% sim - Simulates a sss system using iddata input time series
%
% Syntax:
%       [data,x,tx] = sim(sys,data)
%       [data,x,tx] = sim(sys,data,method)
%       [data,x,tx] = sim(sys,data,method,Ts_sample)
%
% Description:
%       Simulates a sss system using iddata input time series.
%
%       //Note: Requires the System Identification Toolbox.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%       -data:  iddata object containing input time series
%       *Optional Input Arguments:*
%       -method:    time ingetration method
%                   [{'RK4' (continuous) / 'discrete' (discrete)} /
%                   'forwardEuler' / 'backwardEuler' / 'RKDP']
%       -Ts_sample: sampling rate of the state-vector
%
% Output Arguments:
%       -data: iddata object containing output time series
%       -X: matrix of state vectors
%       -tx: time vector for X
%
% Examples:
%       The following code simulates the benchmark "building" (SSS, SISO) 
%       using the backward Euler method
%
%> load building.mat; sys=sss(A,B,C);
%> Ts = 1e-4; %sampling time
%> t = 0:Ts:10; %time vector
%> u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])'; %random gaussian input signal
%> datau = iddata([],u',Ts); %create an iddata object with the input information
%> dataBackward = sim(sys,datau,'backwardEuler'); %simulation with implicit Euler
%> plot(t,dataBackward.y); xlabel('Time [s]'); ylabel('Amplitude');
%
% See Also:
%       iddata/sim, lsim, simForwardEuler, simDiscrete, simBackwardEuler
%
% References:
%       * *[1] Gear (1971)*, Numerical Initial Value Problems in
%       Ordinary Differential Equations.
%       * *[2] Shampine (1994)*, Numerical Solution of Ordinary Differential Equations,
%       Chapman & Hall, New York.
%       * *[3] Shampine and Gordon (1975)*, Computer Solution of Ordinary Differential
%       Equations: the Initial Value Problem, W. H. Freeman, San Francisco.
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
% Authors:      Stefan Jaensch
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if ~any(cellfun(@(x) strcmp('',x),sys.u))
    sys = sys.truncate(':',data.InputName);
else
%     warning('InputName of sys is not specifing. Ignoring InputName definition in iddata')
end

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

if nargout == 1
    y =  lsim(sys,u,Ts,method,Ts_sample);
else
    [y,X,tx] = lsim(sys,u,Ts,method,Ts_sample);
    tx = tx + data.Tstart;
end

data = [data iddata(y',[],Ts,'OutputName',sys.OutputName,'Tstart',data.Tstart)];

