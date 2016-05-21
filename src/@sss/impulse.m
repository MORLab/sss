function  varargout = impulse(varargin)
% IMPULSE - Computes and/or plots the impulse response of a sparse LTI system
%
% Syntax:
%   IMPULSE(sys)
%   IMPULSE(sys,t)
%   IMPULSE(sys1, sys2, ..., t)
%   IMPULSE(sys1,'-r',sys2,'--k',t);
%   [h, t] = IMPULSE(sys,t)
%   [h, t] = IMPULSE(sys,t,opts)
%
% Description:
%       Computes and/or plots the impulse response of a sparse LTI system
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -t:     vector of time values to plot at
%       -Opts:  structure with execution parameters
%			-.odeset:  odeset Settings of ODE solver
%           -.tolOutput: Terminate if norm(y_-yFinal)/norm(yFinal)<tolOutput with yFinal = C*xFinal+D;
%						[1e-3]
%           -.tolState: Terminate if norm(x-xFinal)/norm(xFinal)<tolState with xFinal = -(A\B);
%						[1e-3]
%           -.tf: % return [h, t] instead of tf object as in bult-in case
%                       [0]
%           -.ode: ode solver;
%                       [{'ode45'},'ode113','ode15s','ode23'] 
%
% Outputs:
%       -h, t: vectors containing impulse response and time vector
%
% Examples:
%       The following code computes the impulse response of the benchmark
%       'building' (SSS, SISO) and compares it with the MATLAB built-in function:
%
%> load building.mat
%> sysSparse=sss(A,B,C); %sparse state-space (sss)
%> sys=ss(sysSparse); %full state-space (ss)
%> figure; impulse(sys); hold on; impulse(sysSparse);
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
% Last Change:  10 Nov 2015
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
Def.tf = 0; % return [h, t] instead of tf object as in bult-in case
Def.ode = 'ode45';
Def.nMin = 1000; % impulse responses for models with n<nMin are calculated with the build in Matlab function

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

% final Time
t = [];
tIndex = cellfun(@isfloat,varargin);
if ~isempty(tIndex) && nnz(tIndex)
    t = varargin{tIndex};
    varargin(tIndex)=[];
end

Tfinal = 0;
for i = 1:length(varargin)
    % Set name to input variable name if not specified
    if isprop(varargin{i},'Name')
        if isempty(varargin{i}.Name) % Cascaded if is necessary && does not work
            varargin{i}.Name = inputname(i);
        end
    end
    
    % Convert sss to frequency response data model
    if isa(varargin{i},'sss')
        [TF_] = gettf(varargin{i}, t,Opts);
        varargin{i} = TF_;
        % get length of impulse response
        tMax_ = max(cellfun(@length,TF_.num(:)))*TF_.Ts;
        Tfinal = max(tMax_,Tfinal);
    end
    
end

% Call ss/impulse
if nargout==1 && Opts.tf
    varargout{1} = TF_;
elseif nargout
    [varargout{1},varargout{2},varargout{3},varargout{4}] = impulse(varargin{:},Tfinal);   
    if ~isempty(t)
        if length(t)==1
           t = linspace(0,Tfinal,length(varargout{2})); 
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
else
    if ~isempty(t) && ~isscalar(t)
        % all systems have different Ts (due to ode): plot all systems separately from t(1):sys.Ts:t(end) with impulse-built-in & hold on
        tfindex=zeros(length(varargin),1);
        for i=1:length(varargin)
            if isa(varargin{i},'ss')  ...
                || isa(varargin{i},'tf') || isa(varargin{i},'zpk') ...
                || isa(varargin{i},'frd') || isa(varargin{i},'idtf')...
                || isa(varargin{i},'idpoly') || isa(varargin{i},'idfrd') ...
                || isa(varargin{i},'idss')
                tfindex(i)=1;
            end
        end
        
        for l=1:length(varargin)
            % find indices of the next two systems
            i=find(tfindex,1);
            tfindex(i)=0;
            if ~isempty(i)
                ii=find(tfindex,1);
                if isa(varargin{i},'tf')
                    ti=round(t(1)/varargin{i}.Ts)*varargin{i}.Ts:varargin{i}.Ts:t(end);
                else
                    ti=t;
                end
                if isempty(ii)
                    impulse(varargin{i:end},ti);
                else
                    impulse(varargin{i:ii-1},ti);
                end
            end
            hold on;
        end
        hold off;
    else
        impulse(varargin{:},Tfinal);
    end 
end

function TF = gettf(sys, t, Opts)
Opts.tf = 1;
sys.d = zeros(size(sys.d));
TF = step(sys,t,Opts);



