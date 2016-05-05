function  h = bodeplot(varargin)
% bodeplot - Bode plot of an LTI system
% 
% Syntax: 
%       bodeplot(sys)
%       bodeplot(sys, omega)
%       H = bodeplot(sys)
%       bodemag(sys1, sys2, ..., omega, options)
%
% Description:
%       Bode plot of one or several LTI systems. If an output is defined, a handle to
%       the plot is returned.
%
% Input Arguments:       
%       *Required Input Arguments*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments*
%       -omega:     vector of frequencies or cell with {wmin,wmax}
%       -options:   plot options. see <a href="matlab:help plot">PLOT</a>
%
% Output Arguments:      
%       -h:     plot handle
%
% Examples:
%       Bode plot of the benchmark 'CDplayer' (SSS, MIMO):
%
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> bodemag(sys);
%
% See Also:
%       bode, freqresp, sigma, bodemag
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
% Authors:      Heiko Panzer, Sylvia Cremer, Rudy Eid, 
%               Alessandro Castagnotto
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  06 Apr 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Evaluate options
omegaIndex=0;
for i=1:length(varargin)
    if isa(varargin{i}, 'double')|| isa(varargin{i},'cell')
        omegaIndex = i;
        break
    end
end

if omegaIndex>0
	omega=varargin{omegaIndex};
    varargin(omegaIndex)=[];
else
    omega=[];
end

%% Get frd object
for i=1:length(varargin)
    if isa(varargin{i},'sss')
        varargin{i} = freqresp(varargin{i}, omega, struct('frd',1));
    end
end

%% Call ss/bodeplot
if nargout==1
    h=bodeplot(varargin{:});
else
    bodeplot(varargin{:});
end
