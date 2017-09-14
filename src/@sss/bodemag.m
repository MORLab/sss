function bodemag(varargin)
% BODEMAG - Bode magnitude plot of an LTI system
% 
% Syntax: 
%       BODEMAG(sys)
%       BODEMAG(sys,omega)
%       BODEMAG(sys1,...,sysN,omega)
%       BODEMAG(sys,LineSpec)
%       BODEMAG(sys1,'-r',sys2,'--k', omega);
%
% Description:
%       Bode magnitude plot of one or several LTI systems.
%
% Input Arguments:       
%       *Required Input Arguments*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments*
%       -omega:     vector of frequencies or cell with {wmin,wmax}
%       -LineSpec: String with line style, marker symbol, and color. See <a href="matlab:doc plot">plot</a>
%
% Examples:
%       The following code plots the magnitude of the frequency 
%       response of the benchmark 'CDplayer' (SSS, MIMO):
%
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> bodemag(sys,'-r');
%
% See Also:
%       bode, freqresp, sigma, bodeplot
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
% Authors:      Heiko Panzer, Sylvia Cremer, Rudy Eid, 
%               Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  14 Sep 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Get frd object
varargin = getfrd(varargin{:});

%% Call ss/bodemag
bodemag(varargin{:});
