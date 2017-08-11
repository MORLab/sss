function  varargout = step(varargin)
% STEP - Computes and/or plots the step response of a sparse LTI system
%
% Syntax:
%   STEP(sys)
%   STEP(sys,t)
%   STEP(sys,Tfinal)
%   STEP(sys1, sys2, ..., t)
%   STEP(sys1, sys2, ..., Tfinal)
%   STEP(sys1,'-r',sys2,'--k',t);
%   STEP(sys1,'-r',sys2,'--k',Tfinal)
%   [h, t] = STEP(sys)
%   [h, t] = STEP(sys, t)
%   [h, t] = STEP(sys, Tfinal)
%   [h, t] = STEP(sys, ..., Opts)
%   TF = STEP(sys,...,struct('tf',true))
%   [TF,h,t] = STEP(sys,...,struct('tf',true))
%
% Description:
%       step(sys) plots the step response of the sparse LTI system sys
%
%       [h, t] = step(sys, t) computes the step response of the sparse LTI
%       system sys and returns the vectors h and t with the response and
%       the time series, respectively.
%
%       TF = STEP(sys,struct('tf',true)) returns a discrete time |tf| object 
%       of the FIR-filter with same discrete step response as sys.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -t:     vector of time values to plot at
%       -Tfinal: end time of step response
%       -Opts:  structure with execution parameters
%			-.odeset:  odeset Settings of ODE solver
%           -.tolOutput: Terminate if norm(y_-yFinal)/norm(yFinal)<tolOutput with yFinal = C*xFinal+D;
%						[1e-3 / positive float]
%           -.tolState: Terminate if norm(x-xFinal)/norm(xFinal)<tolState with xFinal = -(A\B);
%						[1e-3 / positive float]
%           -.tf: return tf object
%                       [{0} / 1]
%           -.ode: ode solver;
%                       [{'ode45'} / 'ode113' / 'ode15s' / 'ode23'] 
%           -.tsMin: minimum sample time if no time vector is specified
%                       [{0} / positive float]
%           -.htCell: return ode output as cell with irregularly spaced t
%                       [{0} / 1]
%           -.tLin: uniformly spaced time vector
%                       [{0} / 1]
% 
% Output Arguments:
%       -h, t: vectors containing step response and time vector
%       -TF: discrete time tf object of step response
%
% Examples:
%       The following code plots the step response of the benchmark
%       'building' (SSS, SISO):
%
%> load building.mat; sys=sss(A,B,C);
%> step(sys);
%
% See Also:
%       residue, impulse
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
% Authors:      Heiko Panzer, Sylvia Cremer, Jorge Luiz Moreira Silva
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  14 Jun 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% Make sure the function is used in a correct way before running compts.
nSys = 0;
for iInp = 1:length(varargin)
    if isa(varargin{iInp},'sss') || isa(varargin{iInp},'ss')  ...
            || isa(varargin{iInp},'tf') || isa(varargin{iInp},'zpk') ...
            || isa(varargin{iInp},'frd') || isa(varargin{iInp},'idtf')...
            || isa(varargin{iInp},'idpoly') || isa(varargin{iInp},'idfrd') ...
            || isa(varargin{iInp},'idss')
        nSys = nSys+1;
    end
end
if nSys > 1 && nargout
    error('sss:step:RequiresSingleModelWithOutputArgs',...
        'The "step" command operates on a single model when used with output arguments.');
end

%% Parse inputs and options
Def.odeset = odeset; % Settings of ODE solver
Def.tolOutput = 1e-3; % Terminate if norm(y_-yFinal)/norm(yFinal)<tolOutput with yFinal = C*xFinal+D;
Def.tolState = 1e-3; % Terminate if norm(x-xFinal)/norm(xFinal)<tolState with xFinal = -(A\B);
Def.tf = 0; % return [h, t] instead of tf object as in bult-in case
Def.ode = 'ode45';
Def.nMin = 1000; % impulse responses for models with n<nMin are calculated with the built-in Matlab function
Def.tsMin = 0;
Def.htCell = 0; % return [h,t] cell directly from ode (t is not linearly spaced)
Def.tLin = 0;

% create the options structure
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% store name of the model
name = varargin{1}.Name;

% final Time
t = [];
tIndex = cellfun(@isfloat,varargin);
if ~isempty(tIndex) && nnz(tIndex)
    t = varargin{tIndex};
    varargin(tIndex)=[];
end

if ~isempty(t)
    Tfinal=t(end);
else
    Tfinal=0;
end

for i = 1:length(varargin)
    % Set name to input variable name if not specified
    if isprop(varargin{i},'Name')
        if isempty(varargin{i}.Name) % Cascaded if is necessary && does not work
            varargin{i}.Name = inputname(i);
        end
    end
    
    % Convert sss to frequency response data model
    if isa(varargin{i},'sss') || isa(varargin{i},'ssRed')
        % h,t from ode in cell
        [varargin{i},th] = getht(varargin{i}, t, Opts);
        
        % return [h,t] cell for impulse
        if Opts.htCell
            varargout{1}=varargin{i};
            varargout{2}=th;
            return;
        end

        % get Ts and Tfinal
        Ts=Inf;
        if isempty(t) || isscalar(t)
            for k=1:size(varargin{i},1)
                 for j=1:size(varargin{i},2)
                    Ts=min(min(diff(th{k,j}),Ts));
                    if ~isscalar(t)
                        Tfinal=max(max(th{k,j}),Tfinal);
                    end
                 end
            end
            if Ts<Opts.tsMin
                warning('Ts changed to Opts.tsMin.');
                Ts=Opts.tsMin;
            elseif Ts<1e-5
                warning('Ts is very small. Consider setting Opts.tsMin');
            end
        else
            Ts=min(diff(t));
            Tfinal=t(end);
        end

        if nargout==0
            % create tf object for plotting
            for k=1:size(varargin{i},1)
                 for j=1:size(varargin{i},2)
                     varargin{i}{k,j}=interp1(th{k,j},varargin{i}{k,j},0:Ts:Tfinal);
                     varargin{i}{k,j}(isnan(varargin{i}{k,j}))=0;
                 end
            end
            varargin{i} = cellfun(@(x) [x(1) diff(x)],varargin{i},'UniformOutput',false);
            varargin{i}=filt(varargin{i},1,Ts);
        else
            % [h,t] or tf object output
            if Opts.tf || Opts.tLin
                 % t from 0 to Tfinal with Ts
                 t=0:Ts:Tfinal(end);
            elseif isempty(t) || isscalar(t)
                 % add all time points to a single vector
                 tOut=[];
                 for k=1:size(th,1)
                     for j=1:size(th,2)
                        tOut=[tOut;th{k,j}];
                     end
                 end
                 tOut=sort(tOut);
                 
                 % Tfinal
                 if isscalar(t)
                     tOut=tOut(tOut<Tfinal*ones(size(tOut)));
                 end
                 
                 % keep only elements with abs(u-v)>Ts
                 t=uniquetol(tOut,Ts/tOut(end));
            end
            
            % interpolate h at t
             if nargout ~= 1 || ~Opts.tf
                 h=zeros(length(t),size(varargin{i},1),size(varargin{i},2));
                 for k=1:size(varargin{i},1)
                     for j=1:size(varargin{i},2)
                         h(:,k,j)=interp1(th{k,j},varargin{i}{k,j},t);
                     end
                 end
             end
             
             if nargout==3 && Opts.tf
                 varargin{i} = cellfun(@(x) [x(1) diff(x')],varargin{i},'UniformOutput',false);
                 varargout{1} = filt(varargin{i},1,Ts);
                 varargout{1}.Name = name;
                 varargout{2} = h;
                 varargout{3} = t';
             elseif Opts.tf
                 varargin{i} = cellfun(@(x) [x(1) diff(x')],varargin{i},'UniformOutput',false);
                 varargout{1} = filt(varargin{i},1,Ts);
                 varargout{1}.Name = name;
             else
                 varargout{1}=h;
                 varargout{2}=t;
             end
             return;
        end
    elseif isa(varargin{i},'tf')
         % [h,t] from tf-object -> call built-in
        if isempty(t)
            error('Please specifiy t or Tfinal.');
        end
        [varargout{1},varargout{2},varargout{3},varargout{4}] = step(varargin{:},Tfinal);
        if ~isscalar(t)
            % interpolate [h,0:Ts:Tfinal] at t
            if length(t)==1
               t = linspace(0,t(end),length(varargout{2})); 
            end
            varargout{1} = interp1(varargout{2},varargout{1},t,'spline');
            varargout{2} = t;

            % output uniform with built-in
            if length(size(varargout{1}))==2
                varargout{1}=varargout{1}';
            end
            if size(varargout{2},1)<size(varargout{2},2)
                varargout{2} = varargout{2}';
            end
        end
    end
end

% Plotting with built-in    
if isempty(t)
    t=Tfinal;
end
step(varargin{:},t);

end

function [h,th] = getht(sys, t, Opts)
t = t(:);
h = cell(sys.p,sys.m);
th = cell(sys.p,sys.m);
if sys.n > Opts.nMin
    for i = 1:sys.p
        for j=1:sys.m
            [h{i,j},th{i,j}] = stepLocal(truncate(sys,i,j), t, Opts);
            if size(h{i,j},1)<size(h{i,j},2)
                h{i,j}=h{i,j}';
                th{i,j}=th{i,j}';
            end
        end
    end
else
    for i = 1:sys.p
        for j=1:sys.m
            [h{i,j},th{i,j}] = step(ss(truncate(sys,i,j)), t);
        end
    end
end
end

function [h,te] = stepLocal(sys, t_, Opts)
x0 = zeros(size(sys.x0));
optODE = Opts.odeset;
[A,B,C,D,E,~] = dssdata(sys);
if ~sys.isDescriptor
    odeFun = @(t,x) A*x+B;
else
    % init solveLse
    solveLse(E);
    Opts.reuseLU=true;
    
    % function handle
    odeFun = @(t,x) solveLse(E,A*x+B,Opts);
end
if ~isempty(t_)
    tSim = [0,t_(end)];
else
    xFinal = -(A\B);
    yFinal = C*xFinal+D;
    tSim = [0,decayTime(sys)];
    optODE.OutputFcn = @(t,x,flag) outputFcn(t,x,flag,C,D,xFinal...
        ,yFinal,Opts.tolOutput,Opts.tolState);
end
te = []; h = te;
optODE.Events = @(t,x) eventsFcnT(t,x,C,D);

switch Opts.ode
    case 'ode113'
        [~,~] = ode113(odeFun,tSim,x0,optODE);
    case 'ode15s'
        [~,~] = ode15s(odeFun,tSim,x0,optODE);
    case 'ode23'
        [~,~] = ode23(odeFun,tSim,x0,optODE);
    case 'ode45'
        [~,~] = ode45(odeFun,tSim,x0,optODE);
    otherwise
        error(['The ode solver ' Opts.ode ' is currently not implemented. '  ...
            'Implemented solvers are: ode113, ode15s, ode23, ode45'])
end

    function [value,isterminal,direction]  = eventsFcnT(t,x,C,D)
        value = 1;
        isterminal = 0;
        direction = [];
        y_ = full(C*x+D);
        h = [h, y_'];
        te = [te, t];
    end
end
function status = outputFcn(t,x,flag,C,D,xFinal,yFinal,tolOutput,tolState)
status = 0;
if isempty(x)
    return
end

x = x(:,end);
y_ = full(C*x+D);
if ~isempty(xFinal)
    y_ = full(C*x+D);
    if norm(y_-yFinal)/norm(yFinal)<tolOutput || ...
            norm(x-xFinal)/norm(xFinal)<tolState
        status = 1;
    end
end
end
