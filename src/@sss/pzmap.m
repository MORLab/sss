function [varargout] = pzmap(varargin)
% PZMAP - Pole-zero plot of sparse state-space system
%
% Syntax:
%       PZMAP(sys)
%       PZMAP(sys, opts)
%       [p,z] = PZMAP(sys)
%
% Description:
%       pzmap(sys) creates a pole-zero plot of the continuous- or discrete-time 
%       dynamic system model sys. For MIMO systems, pzmap plots the system 
%       poles and the invariant zeros in one figure. The poles are plotted 
%       as x's and the invariant zeros are plotted as o's.
%
%       [p,z] = pzmap(sys) returns the system poles and invariant zeros in the column 
%       vectors p and z. No plot is drawn on the screen.
%
%//Note: The calculation of the invariant zeros is only defined for systems
%       with the same number of inputs and outputs (m=p). That means that if
%       pzmap is called with a system with m~=p, then z = [ ].
%
% Input Arguments:
%       -sys:      an sss-object containing the LTI system
%       -opts:     plot options. see <a href="matlab:help plot">PLOT</a>
%
% Output Arguments:
%       -p: vector containing poles 
%       -z: vector containing invariant zeros
%
% Examples:
%       Create a random descriptor model (DSSS, SISO) and compare the output
%       of ss/pzmap and sss/pzmap:
%
%> A = randn(500,500); B = randn(500,1); C = randn(1,500); D = zeros(1,1);
%> E = randn(500,500);
%> sys = dss(A,B,C,D,E);
%> sysSss = sss(sys);
%> figure; pzmap(sys);
%> figure; pzmap(sysSss);
%
%       Load the benchmark 'rail_1357' (DSSS, MIMO) and use pzmap:
%
%> load rail_1357.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m),E)
%> figure; pzmap(sys);
%
% See Also:
%       ss/pzmap
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
% Authors:      Heiko Panzer, Sylvia Cremer, Maria Cruz Varona
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  10 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

if isa(varargin{end},'double')
    k=varargin{end};
    varargin=varargin(1:end-1);
else
    k=6;
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
        varargin{i} = zpkData(varargin{i}, k, struct('zpk',true));
    end
end

if nargout
    [varargout{1}, varargout{2}]=pzmap(varargin{:});
else
    pzmap(varargin{:});
end
end

function [varargout] = zpkData(sys,varargin)
% Compute largest poles and zeros or zpk object of an LTI system

%% Parse inputs and options
Def.zpk = false; %return zpk object instead of p and z
Def.typeZ = 'lm'; %eigs type for zeros
Def.typeP = 'lm'; %eigs type for poles

for i=1:length(varargin)
    if isa(varargin{i},'double')
        k=varargin{i};
    elseif isa(varargin{i},'struct')
        Opts=varargin{i};
    end
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

if exist('k','var')
    p=poles(sys,k,struct('type',Opts.typeP));
    z=zeros(sys,k,struct('type',Opts.typeZ));
else
    p=poles(sys,struct('type',Opts.typeP));
    z=zeros(sys,struct('type',Opts.typeZ));
end

if Opts.zpk
    % remove single complex element
    if ~any(isreal(z))
        z(abs(imag(z)-imag(sum(z)))<1e-16)=[];
        p(abs(imag(p)-imag(sum(p)))<1e-16)=[];

        % make sure z and p have the same length
        z=z(1:min(length(z),length(p)));
        p=p(1:length(z));

        if length(z)~=k
            warning(['Because of not complex conjugated values, only ',...
                num2str(length(p)),' zeros and poles are used in the further computation.']);
        end
    end
    
    % k is the first nonzero markov parameter
    if sys.isDae
        error('Does not work for DAE systems yet.');
    else
        k=moments(sys,Inf,2);
        k=k(:,:,2);
    end    
    
    % create zpk object
    varargout{1}=zpk(z,p,k);
else
    varargout{1}=p;
    varargout{2}=z;
end
end
