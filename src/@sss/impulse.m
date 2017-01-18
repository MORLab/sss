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
% Last Change:  18 Jan 2017
% Copyright (c) 2015,2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

    % Make sure the function is used in a correct way before running compts.
    % Check for
    % - output variables allowed only if one system is passed
    % - consistent input/output dimensions throughout models
    nSys = 0;
    m = []; p = [];
    for iInp = 1:length(varargin)
        argClass = class(varargin{iInp});
        if any(strcmp(argClass,{'sss','ss', 'tf', 'zpk', 'frd', 'idtf','idpoly',...
                'idfrd','idss'}))
            % increase system counter
            nSys = nSys+1;
            % check input/output size
            if isempty(m)
                m = size(varargin{iInp}.B,2); p = size(varargin{iInp}.C,1);
            elseif m~= size(varargin{iInp}.B,2) || p~= size(varargin{iInp}.C,1)
                error('sss:impulse:incompatibleInputOutputDimensions',...
                    'When calling sss/impulse with several LTI models, the number of inputs/outputs has to match.');
            end
        end
    end
    if nSys > 1  && nargout
        error('sss:impulse:RequiresSingleModelWithOutputArgs',...
            'The "impulse" command operates on a single model when used with output arguments.');      
    end 

    %% Parse inputs and options
    Def.odeset      = odeset; % Settings of ODE solver
    Def.tolOutput   = 1e-3; % Terminate if norm(y_-yFinal)/norm(yFinal)<tolOutput with yFinal = C*xFinal+D;
    Def.tolState    = 1e-3; % Terminate if norm(x-xFinal)/norm(xFinal)<tolState with xFinal = -(A\B);
    Def.ode         = 'tbd'; %will be determined later depending on the model
    Def.nMin        = 1000; % impulse responses for models with n<nMin are calculated with the built-in Matlab function
    Def.tsMin       = 0; 
    Def.tLin        = 0; % linearly spaced t

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
    g = cell(1,nSys); tg = cell(1,nSys); % Preallocation
    
    for iArgin = 1:length(varargin) % Loop through all systems
        
        % Set name to input variable name if not specified
        if isprop(varargin{iArgin},'Name')
            if isempty(varargin{iArgin}.Name) % Cascaded if is necessary && does not work
                varargin{iArgin}.Name = inputname(iArgin);
            end
        end

        % Convert sss to frequency response data model
        if isa(varargin{iArgin},'sss')        
            % get g,t from ode
            if ~isempty(t)
                [g{iArgin},tg{iArgin}] = getgt(varargin{iArgin}, t(end),Opts);
            else
                [g{iArgin},tg{iArgin}] = getgt(varargin{iArgin}, [], Opts);
            end
        elseif isa(varargin{iArgin},'tf') || isa(varargin{iArgin},'ss')
            % get g,t from built-in impulse
            [g{iArgin},tg{iArgin}] = impulse(varargin{:},Tfinal);
        end
    end
    
    %% Plotting or returning outputs
    if nargout
        varargout{1} = g{1};
        varargout{2} = tg{1};
    else
        %== Initialize plot as ltiplot (resppack.timeplot)
        ax = gca;
        h = ltiplot(ax,'impulse',varargin{1}.InputName,varargin{1}.OutputName,[],cstprefs.tbxprefs);
        set(h,  'Visible','on')
        
        %== Get axes handles to subplots
        ah = get(gcf,'Children');
        ah = ah(strcmp('axes',get(ah,'type')));
        if length(ah)~= m*p
            error('sss:impulse:wrongNumberOfSubplots','Something went wrong in generating or detecting the subplots for the impulse response');
        end
        set(ah, 'XLimMode','auto','YLimMode','auto');

        %== plot computed data       
        for kSys = 1:nSys
            cnt = 0;
            for iIn = m:-1:1
                for jOut = p:-1:1
                    cnt = cnt+1;
                    line(tg{kSys},g{kSys}(:,iIn,jOut),...
                        'Parent',ah(cnt),'DisplayName',varargin{kSys}.Name,...
                        'Color',pickColor(kSys));
                end
            end
        end
    end
end

%% == AUXILIARY

function [g,te] = getgt(sys, t_, Opts)
    % In:
    % sys: system to be simulated
    % t_:  initial/final time
    % Opts
    
    % Out:
    % te: l-dimensional time vector
    % g: lxpxm array, where p is sys.p and m is sys.m
    
    % Define ODE options
    optODE  = Opts.odeset;
    optODE  = odeset(optODE,'Vectorized','on','Jacobian',sys.A);
    
    % Define ODE fun
    if ~sys.isDescriptor
        odeFun = @(t,x) sys.A*x;
    else
        % init solveLse
        solveLse(sys.E);
        Opts.reuseLU=true;

        % function handle
        odeFun = @(t,x) solveLse(sys.E,sys.A*x,Opts); %more efficient than built-in for sparse models
    end
    
    if ~isempty(t_)
        tSim = [0,t_(end)];
    else
        tSim = [0,decayTime(sys)];
    end
    
    xFinal = zeros(sys.n,1);
    optODE.OutputFcn = @(t,x,flag) outputFcn(t,x,flag,sys.C,sys.D,xFinal,...
        Opts.tolOutput,Opts.tolState);    
    
    te = []; g = [];

    if strcmp(Opts.ode,'tbd')
        solver = chooseSolver(sys);
    else
        solver = Opts.ode;
    end

    for iIn = 1:sys.m
        x0 = sys.B(:,iIn); %take respective column of B as initial condition
        switch solver
            case 'ode113', ode113(odeFun,tSim,x0,optODE);
            case 'ode15s', ode15s(odeFun,tSim,x0,optODE);
            case 'ode23',  ode23(odeFun,tSim,x0,optODE);
            case 'ode45',  ode45(odeFun,tSim,x0,optODE);
            otherwise
                error(['The ode solver ' Opts.ode ' is currently not implemented. '  ...
                    'Implemented solvers are: ode113, ode15s, ode23, ode45'])
        end
        te = cat(3,te,tCurr);
        g  = cat(3,g, gCurr);
    end

%     function [value,isterminal,direction]  = eventsFcnT(t,x,C,D)
%         value = 1;
%         isterminal = 0;
%         direction = [];
%     end
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
        ode = 'ode15s'; 
%         ode = 'ode45'; 
    end

    function status = outputFcn(t,x,flag,C,D,xFinal,tolOutput,tolState)
        status = 0; %initialized as not converged
        if isempty(x)
            return
        end
        
        % Generate the current output
        y_  = full(C*x);
        
        %store only the solution (instead of the whole x vector
        
        % tCurr: a lx1 array of time instances
        % gCurr: a lxp array of outputs
        if strcmp(flag,'init')
            gCurr  = y_.';
            tCurr  = t(1);
        else
            gCurr = [gCurr; y_.']; 
            tCurr = [tCurr; t.'];
        end
    
        % Check for convergence to the final solution
        if ~isempty(xFinal)
            yFinal = C*xFinal;
            if norm(y_(:,end)-yFinal)/norm(yFinal)<tolOutput || ...
                    norm(x(:,end)-xFinal)/norm(xFinal)<tolState
                status = 1;
            end
        end
    end

end

function c = pickColor(counter)
    % Trivial function to loop through plot colors
    
    %   Default colors
    C= [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
    
    nC = size(C,1); %number of colors
    
    % the actual piece of magic
    currC = counter - nC*floor((counter-1)/nC);    
    c = C(currC,:);
end






