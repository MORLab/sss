function  [varargout] = bode(varargin)
% BODE Plots the bode diagram of an LTI system
%
%
% Syntax:
%   BODE(sys)
%   BODE(sys,omega)
%   BODE(sys1, sys2, ..., omega)
%   BODE(sys1,'-r',sys2,'--k',w);
%   [mag, phase, omega] = BODE(sys)
%   frdData = BODE(sys,...,struct('frd',1))
%
% Description:
%       This function computes the bode plot of one or several LTI systems
%       given as sss objects. If no ouput is specified, then a plot is
%       generated.
%
%       If the frequency range is not specified, the function will
%       determine it. It is also possible to pass several systems or
%       plotting options, just like in MATLAB's built-in version.
%
%       It the function is called with only one ouput and the option 'frd'
%       is specified as last input variable, than an frd object is
%       returned.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:   sss-object containing the LTI system
%       -omega: a vector of frequencies
%       *Optional Input Arguments:*
%       -Opts:  structure with execution parameters
%			-.frd:  return frd object;
%						[{0} / 1]%
% Output Arguments:
%       - mag/phase: magnitude and phase response
%       - omega:     frequencies corresponding to the data
%       - frdData:   a frd object with the frequency response data
%
% Examples:
%		This code loads a benchmark model included in the toolbox
%		and plots its bode diagram using the sparse state-space class:
%
%> load building; 
%> sys = sss(A,B,C);
%> bode(sys);
%
% See Also:
%   freqresp, sigma
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
% Authors:      Heiko Panzer, Stefan Jaensch, Sylvia Cremer,
%               Rudy Eid, Alessandro Castagnotto, Lisa Jeschek
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  21 Mar 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% Parse inputs and options
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
else
    Opts = struct();
end

Def.frd = 0; %return magnitude instead of frd object as in bult-in case 

% create the options structure
if ~exist('Opts','var') || isempty(fieldnames(Opts))
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% Make sure the function is used in a correct way before running compts.
nSys = 0;
for iInp = 1:length(varargin)
    if isa(varargin{iInp},'sss')
        nSys = nSys+1;
    end   
end
if nSys > 1 && nargout
    error('sss:bode:RequiresSingleModelWithOutputArgs',...
        'The "bode" command operates on a single model when used with output arguments.');
end

% Frequency vector
omega = [];
omegaIndex = cellfun(@isfloat,varargin);
if ~isempty(omegaIndex) && nnz(omegaIndex)
    omega = varargin{omegaIndex};
    varargin(omegaIndex)=[];
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
        varargin{i} = getfrd(varargin{i}, omega);
    end
end

% Call ss/bode
if nargout == 1 && Opts.frd
    varargout{1} = varargin{1};
elseif nargout
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}] = bode(varargin{:});
else
    bode(varargin{:});
end
end

function resp = getfrd(sys, omega)

if not(exist('omega','var')) || isempty(omega)
    [m, w] = freqresp(sys);
else
    [m, w] = freqresp(sys,omega);
end

%  remove frequencies at infinity to create frd object
k = find(isinf(w));
w(k) = []; m(:,:,k) = [];

resp = frd(m,w,sys.Ts,...
    'InputName',sys.InputName,'OutputName',sys.OutputName,...
    'Name',sys.Name);
end