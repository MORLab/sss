function  varargout = impulse(varargin)
% IMPULSE - Computes and/or plots the impulse response of a sparse LTI system
%
% Syntax:
%   IMPULSE(sys)
%   IMPULSE(sys,t)
%   IMPULSE(sys,Tfinal)
%   IMPULSE(sys1, sys2, ..., t)
%   IMPULSE(sys1, sys2, ..., Tfinal)
%   IMPULSE(sys1,'-r',sys2,'--k',t)
%   IMPULSE(sys1,'-r',sys2,'--k',Tfinal)
%   [y, t] = IMPULSE(sys,t)
%   [y, t] = IMPULSE(sys,Tfinal)
%   [y, t] = IMPULSE(sys,...,Opts)
%   TF     = IMPULSE(sys,...,struct('tf',true))
%   [TF,y,t] = IMPULSE(sys,...,struct('tf',true))
%
% Description:
%       Computes and/or plots the impulse response of a sparse LTI system
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -t:     vector of time values to plot at
%       -Tfinal: end time of impulse response
%       -Opts:  structure with execution parameters
%			-.odeset:  settings of ODE solver
%           -.tolOutput: Terminate if norm(y_-yFinal)/norm(yFinal)<tolOutput with yFinal = C*xFinal+D;
%						[1e-3 / positive float]
%           -.tolState: Terminate if norm(x-xFinal)/norm(xFinal)<tolState with xFinal = -(A\B);
%						[1e-3 / positive float]
%           -.tf: return tf object
%                       [{0} / 1]
%           -.ode: ode solver;
%                       [{'tbd} / 'ode45' / 'ode113' / 'ode15s' / 'ode23'] 
%           -.tsMin: minimum sample time if no time vector is specified
%                       [{0} / positive float]
%           -.tLin: uniformly spaced time vector
%                       [{0} / 1]
%
% Outputs:
%       -y, t: vectors containing impulse response and time vector
%       -TF: discrete time |tf| object of step response
%
% Examples:
%       The following code computes the impulse response of the benchmark
%       'building' (SSS, SISO) and compares it with the MATLAB built-in function:
%
%> load beam.mat
%> sys=sss(A,B,C);    %sparse state-space (sss)
%> figure; impulse(sys); 
%> hold on; impulse(ss(sys)); %compare to built-in
%> legend('ss/impulse','sss/impulse');
%
% See Also:
%       residue, ss/impulse, step
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
% Authors:      Heiko Panzer, Sylvia Cremer, Jorge Luiz Moreira Silva
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
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
        error('sss:impulse:RequiresSingleModelWithOutputArgs',...
            'The "impulse" command operates on a single model when used with output arguments.');
    end

    %% Parse inputs and options
    Def.odeset = odeset; % Settings of ODE solver
    Def.tolOutput = 1e-3; % Terminate if norm(y_-yFinal)/norm(yFinal)<tolOutput with yFinal = C*xFinal+D;
    Def.tolState = 1e-3; % Terminate if norm(x-xFinal)/norm(xFinal)<tolState with xFinal = -(A\B);
    Def.tf = 0; % return [y, t] instead of tf object as in bult-in case
    Def.ode = 'tbd'; %will be determined later depending on the model
    Def.nMin = 1000; % impulse responses for models with n<nMin are calculated with the built-in Matlab function
    Def.tsMin = 0; 
    Def.tLin =0; % linearly spaced t

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


    %% store name of the model
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

    %%  Computation
    for i = 1:length(varargin)
        % Set name to input variable name if not specified
        if isprop(varargin{i},'Name')
            if isempty(varargin{i}.Name) % Cascaded if is necessary && does not work
                varargin{i}.Name = inputname(i);
            end
        end

        % Convert sss to frequency response data model
        if isa(varargin{i},'sss')
            % get g,t from ode in cell
            if ~isempty(t)
                [varargin{i},tg] = getgt(varargin{i}, t(end),Opts);
            else
                [varargin{i},tg] = getgt(varargin{i}, [], Opts);
            end
keyboard %TODO: continue debugging from here

    %         % compute impulse response [g,t]
    %         for k=1:size(varargin{i},1)
    %             for j=1:size(varargin{i},2)
    %                 varargin{i}{k,j}=gradient(varargin{i}{k,j},tg{k,j});
    %                 tg{k,j}(isnan(varargin{i}{k,j}))=[];
    %                 varargin{i}{k,j}(isnan(varargin{i}{k,j}))=[];
    %             end
    %         end

            % get Ts and Tfinal
            Ts=Inf;
            if isempty(t) || isscalar(t)
                for k=1:size(varargin{i},1)
                     for j=1:size(varargin{i},2)
                        Ts=min(min(diff(tg{k,j}),Ts));
                        if ~isscalar(t)
                            Tfinal=max(max(tg{k,j}),Tfinal);
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
                % compute tf object for plotting
                for k=1:size(varargin{i},1)
                     for j=1:size(varargin{i},2)
                         varargin{i}{k,j}=Ts*interp1(tg{k,j},varargin{i}{k,j},0:Ts:Tfinal);
                         varargin{i}{k,j}(isnan(varargin{i}{k,j}))=0;
                     end
                end
                varargin{i}=filt(varargin{i},1,Ts);
            else
                 % [g,t] or tf object output
                 if Opts.tf || Opts.tLin
                     % t from 0 to Tfinal with Ts
                     t=0:Ts:Tfinal(end);
                 elseif isempty(t) || isscalar(t)
                     % add all time points to a single vector
                     tOut=[];
                     for k=1:size(tg,1)
                         for j=1:size(tg,2)
                            tOut=[tOut;tg{k,j}];
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

                 % interpolate g at t
                 if nargout ~= 1 || ~Opts.tf
                     g=zeros(length(t),size(varargin{1},1),size(varargin{1},2));
                     for k=1:size(varargin{1},1)
                         for j=1:size(varargin{1},2)
                             g(:,k,j)=interp1(tg{k,j},varargin{1}{k,j},t);
                         end
                     end
                 end

                 if Opts.tf
                    % compute tf object
                    for k=1:size(varargin{i},1)
                         for j=1:size(varargin{i},2)
                             varargin{i}{k,j}=Ts*interp1(tg{k,j},varargin{i}{k,j},t);
                             varargin{i}{k,j}(isnan(varargin{i}{k,j}))=0;
                         end
                    end
                    TF=filt(varargin{i},1,Ts);
                 end

                 if nargout==3 && Opts.tf
                     varargout{1} = TF;
                     varargout{1}.Name = name;
                     varargout{2} = g;
                     varargout{3} = t';
                 elseif Opts.tf
                     varargout{1} = TF;
                     varargout{1}.Name = name;
                 else
                     varargout{1}=g;
                     varargout{2}=t;
                 end
                 return;
            end
        elseif isa(varargin{i},'tf') && nargout
            % [g,t] from tf-object -> call built-in
            if isempty(t)
                error('Please specifiy t or Tfinal.');
            end
            [varargout{1},varargout{2},varargout{3},varargout{4}] = impulse(varargin{:},Tfinal);
            if ~isscalar(t)
                % interpolate [g,0:Ts:Tfinal] at t
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

    %% Plotting with built-in
    if isempty(t)
        t=Tfinal;
    end
    impulse(varargin{:},t);

end

function [g,tg] = getgt(sys, te, Opts)
    Opts.tf = false;
    Opts.htCell = true;
     
    sys.d = zeros(size(sys.d));
    
    t = te;
    g = cell(sys.p,sys.m);
    tg = cell(sys.p,sys.m);
    for i = 1:sys.p
        for j=1:sys.m
            [g{i,j},tg{i,j}] = impulseLocal(truncate(sys,i,j), t, Opts);
            if size(g{i,j},1)<size(g{i,j},2)
                g{i,j}=g{i,j}';
                tg{i,j}=tg{i,j}';
            end
        end
    end
end
    
function [g,te] = impulseLocal(sys, t_, Opts)
    
    x0 = sys.B;
    optODE = Opts.odeset;
    [A,B,C,D,E,~] = dssdata(sys);
    if ~sys.isDescriptor
        odeFun = @(t,x) A*x;
    else
        % init solveLse
        solveLse(E);
        Opts.reuseLU=true;

        % function handle
        odeFun = @(t,x) solveLse(E,A*x,Opts);
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
    te = []; g = te;
    optODE.Events = @(t,x) eventsFcnT(t,x,C,D);

    if strcmp(Opts.ode,'tbd')
        solver = chooseSolver(sys);
    else
        solver = Opts.ode;
    end

    switch solver
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
        g = [g, y_'];
        te = [te, t];
    end
    function ode = chooseSolver(sys)
        % choose between:
        % ode23, ode45, ode113, ode15s, ode23s, ode23t, and ode23tb
        %
        % Characteristics: 
        %{
            ode45  | Nonstiff | Medium | Most of the time. This should be the first solver you try.
            ode23  | Nonstiff | Low    | For problems with crude error tolerances or for solving moderately stiff problems.
            ode113 | Nonstiff | Low to high | For problems with stringent error tolerances or for solving computationally intensive problems.
            ode15s | Stiff    | Low to medium | If ode45 is slow because the problem is stiff.
            ode23s | Stiff    | Low     | If using crude error tolerances to solve stiff systems and the mass matrix is constant.
            ode23t | Moderately Stiff | Low | For moderately stiff problems if you need a solution without numerical damping.
            ode23tb| Stiff | Low | If using crude error tolerances to solve stiff systems.
        %}
        
        %hard coded. Try to find out if the problem is stiff
%         ode = 'ode15s'; 
        ode = 'ode45'; 
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




