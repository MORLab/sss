function sysfrd = frd(varargin)
% FRD - Convert to frequency-response data model
% 
% Syntax:
%       sysfrd = FRD(sys)
%       sysfrd = FRD(sys, omega)
%       sysfrd = FRD(sys, Opts)
%       sysfrd = FRD(sys, omega, Opts)
%
% Description:
%       Evaluates complex transfer function of LTI systems and return an
%       frd object.
%
%       As opposed to the built-in function, |frd(sys)| can be called without
%       specifying a frequency vector, which will then be determined
%       automatically.
%
% Inputs:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -omega: vector of frequencies or cell with {wmin,wmax}
%       -Opts:  structure with execution parameters
%           -.maxPoints: Maximum number of refinement points
%                       [{1500} / positive integer]
%           -.lse:  solve linear system of equations
%                       [{'sparse'} / 'full' /'gauss' /'hess' / 'iterative']
%       
% Outputs:      
%       -sysfrd: a MATLAB built-in frdobject with the frequency-response
%
% Examples:
%       The following code computes the frequency response of the benchmark
%       'building' and returns a respective frd object, which is then 
%       plotted using |bodemag|.
%
%> sys = sss('building.mat');
%> sysfrd=frd(sys);
%> bodemag(sysfrd);
%
% See Also:
%       freqresp, bode, sigma, bodemag, bodeplot
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
% Authors:      Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  14 Sep 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Input parsing
Def.maxPoints   = 1500; % maximum number of refinement points
Def.lse         = 'sparse'; % solveLse

narginchk(1,3)

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

% define sys and omega
sys = varargin{1};
if length(varargin)>1, omega = varargin{2}; else omega = [];end
       
%% Run computations

[G, omega] = freqresp(sys,omega,Opts);

%  remove frequencies at infinity to create frd object
k = find(isinf(omega));
omega(k) = []; G(:,:,k) = [];

sysfrd = frd(G,omega,sys.Ts,...
            'InputName',sys.InputName,'OutputName',sys.OutputName,...
            'Name',sys.Name);

