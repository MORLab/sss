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
%   [h, t] = IMPULSE(sys,t)
%   [h, t] = IMPULSE(sys,Tfinal)
%   [h, t] = IMPULSE(sys,...,Opts)
%   TF     = IMPULSE(sys,...,struct('tf',true))
%   [TF,h,t] = IMPULSE(sys,...,struct('tf',true))
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
%
% Outputs:
%       -h, t: vectors containing impulse response and time vector
%       -tf: discrete time tf object of step response
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
Def.tf = 0; % return [h, t] instead of tf object as in bult-in case
Def.ode = 'ode45';
Def.nMin = 1000; % impulse responses for models with n<nMin are calculated with the build in Matlab function
Def.tsMin = 0; 

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
Tfinal=[];
tIndex = cellfun(@isfloat,varargin);
if ~isempty(tIndex) && nnz(tIndex)
    t = varargin{tIndex};
    varargin(tIndex)=[];
end

for i = 1:length(varargin)
    % Set name to input variable name if not specified
    if isprop(varargin{i},'Name')
        if isempty(varargin{i}.Name) % Cascaded if is necessary && does not work
            varargin{i}.Name = inputname(i);
        end
    end
    
    % Convert sss to frequency response data model
    if isa(varargin{i},'sss')        
        if isempty(t) || isscalar(t)
            [varargin{i},Tmax] = gettf(varargin{i}, t, Opts);
            Tfinal=max([Tmax,Tfinal]);
        else
            Ts=min(diff(t));
            [h,th] = getht(varargin{i}, t(end),Opts);
            g=cell(varargin{i}.p,varargin{i}.m);
            for k=1:varargin{i}.p
                for j=1:varargin{i}.m
                    if size(h{k,j},1)>size(h{k,j},2)
                        h{k,j}=h{k,j}';
                        th{k,j}=th{k,j}';
                    end
                    g{k,j}=Ts*interp1(th{k,j}(1:end-1)+diff(th{k,j})/2,diff(h{k,j})./diff(th{k,j}),0:Ts:t(end));
                    g{k,j}(isnan(g{k,j}))=0;
                end
            end
            varargin{i}=filt(g,1,Ts);
        end
    end
    
end

if isempty(t)
    t=Tfinal;
end

% Call ss/impulse
if nargout==1 && Opts.tf
    varargout{1} = varargin{1};
elseif nargout==3 && Opts.tf
    varargout{1} = varargin{1};
    Ts = varargin{1}(1,1).Ts;
    n = length(varargin{1}(1,1).num{1})-1;
    varargout{2} = cell(size(varargin{1}));
    for i=1:size(varargin{1},1)
        for j=1:size(varargin{1},2) 
            varargout{2}{i,j} = varargin{1}(i,j).num{1,1} / Ts;
        end
    end
    varargout{3} = 0:Ts:n*Ts;
    if size(varargout{2},2)==1 && size(varargout{2},1)==1
       varargout{2} = varargout{2}{1,1};
    end
elseif nargout
    [varargout{1},varargout{2},varargout{3},varargout{4}] = impulse(varargin{:},t(end));   
    if ~isempty(t)
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
else
    impulse(varargin{:},t);
end

function [TF,tMax] = gettf(sys, t, Opts)
Opts.tf = true;
Opts.htOde = false;
sys.d = zeros(size(sys.d));
[TF,h,~] = step(sys,t,Opts);
tMax = max(cellfun(@length,TF.num(:)))*TF.Ts;

function [h,t] = getht(sys, te, Opts)
Opts.tf = false;
Opts.htOde = true;
sys.d = zeros(size(sys.d));
[h,t] = step(sys,te,Opts);



