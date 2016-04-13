function  varargout = impulse(varargin)
% IMPULSE - Computes and/or plots the impulse response of a sparse LTI system
%
% Syntax:
%       [h, t] = impulse(sys,t)
%       [h, t] = impulse(sys,t,opts)
%
% Description:
%       Computes and/or plots the impulse response of a sparse LTI system
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -t:     vector of time values to plot at
%       -opts:  plot options. see <a href="matlab:help plot">PLOT</a>
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

%% Parse inputs and options
Def.odeset = odeset;
Def.tolOutput = 1e-3; % Terminate if norm(y_-yFinal)/norm(yFinal)<tolOutput with yFinal = C*xFinal+D;
Def.tolState = 1e-3; % Terminate if norm(x-xFinal)/norm(xFinal)<tolState with xFinal = -(A\B);
Def.tf = 0; %return [h, t] instead of tf object as in bult-in case
Def.ode = 'ode45';
Def.nMin = 1000;

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
else
    impulse(varargin{:},Tfinal);
end

function TF = gettf(sys, t, Opts)
Opts.tf = 1;
TF = step(sys,t,Opts);

TF.Ts


