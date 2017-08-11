function [varargout] = pzmap(varargin)
% PZMAP - Pole-zero plot of sparse state-space system
%
% Syntax:
%       PZMAP(sys)
%       PZMAP(sys, k)
%       PZMAP(sys, ..., Opts)
%       [p,z] = PZMAP(sys, k, Opts)
%
% Description:
%       pzmap(sys) creates a pole-zero plot of the continuous- or discrete-time 
%       dynamic system model sys containing only the 6 poles and zeros with
%       largest magnitude. This number can be changed with input k. The type
%       of the computed poles and zeros can be specified with the options 'typeP'
%       and 'typeZ'. For MIMO systems, pzmap plots the system poles and the 
%       invariant zeros in one figure. The poles are plotted as x's and the 
%       invariant zeros are plotted as o's.
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
%       -k:        number of computed poles and zeros
%       -Opts:     structure with execution parameters
%			-.typeP:    type of poles
%						[{'lm'} / 'sm' / 'la' / 'sa']
%			-.typeZ:    type of zeros
%						[{'lm'} / 'sm' / 'la' / 'sa']
%
% Output Arguments:
%       -p: vector containing poles 
%       -z: vector containing invariant zeros
%
% Examples:
%       Load the benchmark 'building' (SSS, SISO) and use pzmap to plot the
%       first poles and zeros with largest magnitude. Compare the result to
%       the built-in function that uses dense (!) operations to compute the
%       whole spectrum.
%
%> load building.mat, sys = sss(A,B,C);
%> figure; pzmap(ss(sys));hold on; pzmap(sys);
%
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
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Sylvia Cremer, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  10 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% Parse inputs and options
Def.typeZ = 'lm'; %eigs type for zeros
Def.typeP = 'lm'; %eigs type for poles

if isa(varargin{end},'struct')
    Opts=varargin{end};
    varargin=varargin(1:end-1);
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

k=6;
i=2;
while(i<=length(varargin))
    if isa(varargin{i},'double')
        k=varargin{i};
        varargin(i)=[];
    end
    i=i+1;
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
        if varargin{i}.isDae
            error('pzmap does not work with DAE systems yet.');
        end
        varargin{i} = zpk(varargin{i}, k, Opts.typeP, Opts.typeZ);
    end
end

if nargout
    [varargout{1}, varargout{2}]=pzmap(varargin{:});
else
    pzmap(varargin{:});
end
end


