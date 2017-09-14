function [varargout] = bode(varargin)
% BODE - Plots the bode diagram of an LTI system
%
% Syntax:
%   BODE(sys)
%   BODE(sys, omega)
%   BODE(sys1, sys2, ..., sysN)
%   BODE(sys1, sys2, ..., sysN, omega)
%   BODE(sys1,'-r',sys2,'--k');
%   BODE(sys1,'-r',sys2,'--k',omega);
%   [mag, phase, omega] = BODE(sys)
%   [mag, phase, omega] = BODE(sys, omega)
%
% Description:
%       This function computes the bode plot of one or several LTI systems
%       given as sss objects. If no ouput is specified, then a plot is
%       generated. 
%
%       The first model passed to BODE needs to be an sss object, any other 
%       may be a built-in object of classes ss and frd.
%
%       If the frequency range is not specified, the function will
%       determine it. It is also possible to pass several systems or
%       plotting options, just like in MATLAB's built-in version.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:   sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -omega: vector of frequencies or cell with {wmin,wmax}
%       -LineSpec: String with line style, marker symbol, and color. See <a href="matlab:doc plot">plot</a>
%
% Output Arguments:
%       - mag/phase: magnitude and phase response
%       - omega:     frequencies corresponding to the data
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
%   freqresp, sigma, bodeplot, bodemag
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
% Authors:      Heiko Panzer, Stefan Jaensch, Sylvia Cremer, Rudy Eid,
%               Alessandro Castagnotto, Lisa Jeschek, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  14 Sep 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% Evaluate options
sssIndex=0;
for i=1:length(varargin)
    if isa(varargin{i},'sss') || isa(varargin{i},'ssRed')
        sssIndex = i;
    end
    if sssIndex > 1, break; end
end

if nargout
    if sssIndex > 1 
        error('sss:bode:RequiresSingleModelWithOutputArgs',...
            'The "bode" command operates on a single model when used with output arguments.');
    end
end

%% Get frd object
varargin = getfrd(varargin{:});

%% Call ss/bode or ss/bodeplot
if nargout
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}] = bode(varargin{:});
else
    bodeplot(varargin{:});
end
end