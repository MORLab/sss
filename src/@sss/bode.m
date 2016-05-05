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
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:   sss-object containing the LTI system
%       -omega: vector of frequencies or cell with {wmin,wmax}
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
% Last Change:  06 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

% Make sure the function is used in a correct way before running compts.
for iInp = 1:length(varargin)
    if isa(varargin{iInp},'double') || isa(varargin{iInp},'cell')
        omegaIndex=iInp;
    end
end

% Frequency vector
if exist('omegaIndex','var') && nnz(omegaIndex)
    omega = varargin{omegaIndex};
    varargin(omegaIndex)=[];
end

if length(varargin) > 1 && nargout
    error('sss:bode:RequiresSingleModelWithOutputArgs',...
        'The "bode" command operates on a single model when used with output arguments.');
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
        if not(exist('omega','var')) || isempty(omega)
            varargin{i} = freqresp(varargin{i},struct('frd',1));
        else
            varargin{i} = freqresp(varargin{i},omega, struct('frd',1));
        end
    end
end

% Call ss/bode
if nargout
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}] = bode(varargin{:});
else
    bodeplot(varargin{:});
end
end