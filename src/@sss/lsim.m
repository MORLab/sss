function [y,tx,x] = lsim(sys,u,t,x0,Opts)
% LSIM - Simulate time response of sparse dynamic system to arbitrary inputs
%
% Syntax:
%       y        = LSIM(sys,u,t)
%       [y,x,tx] = LSIM(sys,u,t,x0)
%       [y,x,tx] = LSIM(sys,u,t,x0,Opts)
%
% Description:
%       Simulates the response in time of sparse LTI systems to arbitrary inputs.
%
%       The initial state |x0| can be passed either as an input argument or
%       by defining the respective property |sys.x0| in the sss object.
%
%       The time variable |t| can be either a vector of the form 0:ts:tEnd,
%       just as in the built-in case, or directly a scalar with the
%       sampling time ts.
%
%       Optional parameters such as a different sampling time for the state
%       vector or the simulation method can be defined through the option
%       structure |Opts|.
%       
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%       -u:     Input signal as vector/matrix (Nsamples x Ninputs)
%       -t:     Time vector for u (alternatively: sampling time)
%
%       *Optional Input Arguments:*
%       -Opts:	structure with execution parameters
%			-.method:    time ingetration method;
%                       [{'RK4' (continuous) / 'discrete' (discrete)} / 'forwardEuler' / 'backwardEuler' / 'RKDP']
%           -.Ts_sample: sampling rate of the state-vector;
%                       [{Ts}, float]
%
% Output Arguments:
%       -y:     Matrix with the response of the system  (Nsamples x Noutputs)
%       -tx:    time vector for X                       (Nsamples x 1)
%       -x:     Matrix of state vectors                 (Nsamples x Nstates)
%
% Examples:
%		Simulate the response of the building model to a ramp input.
%
%> sys      = sss('building');
%> u        =  [0:.1:100].';
%> ts       = 0.01; %sampling time
%> [y,t,x]  = lsim(sys,u,ts);
%
%       Check the results, e.g. by plotting the output and the second state
%       variable over time.
%> figure; plot(t,y)
%
%> figure; plot(t,x(:,2))
%
% See Also:
%       sss/sim sss/lsim
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
% Authors:      Michael Leipold, Stefan Jaensch, Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  04 Aug 2018
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


%% Input parsing
narginchk(3,5);

% Define default execution parameters
if sys.Ts == 0 %continuous time model
    Def.method = 'RK4';
else
    % discrete time model
    Def.method = 'discrete';
end

% Parse time vector
if length(t) ~= size(u,1)
    if isscalar(t) %sampling time was passed
        Ts  = t;
%         t   = 0:Ts:Ts*(size(u,1)-1);
    else
        error('sss:lsim:timeVector','The time vector t must have the same length as the rows of u') 
    end
else %time vector passed.
    %Check t for equidistant steps and get sampling time
    Ts = mean(diff(t)); tol = 1e-10; %numerical rounding errors
    if Ts - max(diff(t)) > tol || Ts - min(diff(t)) > tol
        error('sss:lsim:timeVector','The time vector t must be of the form 0:ts:tEnd') 
    end
end
Def.Ts_sample = Ts;

% create the options structure
if ~exist('Opts','var') || isempty(fieldnames(Opts))
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% Define initial state vector
if nargin < 4 || isempty(x0)
    % take from sys object
    x0 = sys.x0;
end

%% Perform simulation

[A,B,C,D,E] = dssdata(sys);

simInputs = {A,B,C,D,E,u,x0,Ts,Opts.Ts_sample,sys.isDescriptor};

switch Opts.method
    case 'forwardEuler'
        if nargout == 1
            y        = simForwardEuler(simInputs{:});
        else
            [y,x,tx] = simForwardEuler(simInputs{:});
        end
    case 'backwardEuler'
        if nargout == 1
            y           = simBackwardEuler(simInputs{:});            
        else
            [y,x,tx]    = simBackwardEuler(simInputs{:});
        end
    case 'RK4'
        if nargout == 1
            y           = simRK4(simInputs{:});
        else
            [y,x,tx]    = simRK4(simInputs{:});
        end
    case 'discrete'
        if nargout == 1
            y           = simDiscrete(simInputs{:});
        else
            [y,x,tx]    = simDiscrete(simInputs{:});
        end
end

%% Transpose y to have the same size as built-in lsim
y   = y.';
if nargout > 1
    tx  = tx.';
    x   = x.';
end

end

