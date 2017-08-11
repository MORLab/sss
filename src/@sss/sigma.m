function  [s, omega] = sigma(varargin)
% sigma - Plots the singular values of the frequency response of a sparse LTI system
% 
% Syntax: 
%       sigma(sys)
%       sigma(sys, omega)
%       [s, omega] = sigma(sys)
%       sigma(sys1,...,sysN)
%       sigma(sys1,...,sysN, omega)
%       sigma(sys1,LineSpec,...,sysN,omega)
%
% Description:
%       Computes and plots the singular values of the frequency response of one or 
%       several sparse LTI systems. The frequency vector is automatically
%       chosen by calling |freqresp|, if not provided. 
%
%       The function allows the usage of several |sss| and |ss| objects,
%       including optional LineSpecs.
%
% Input Arguments:       
%       *Required Input Arguments*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments*
%       -omega:     vector of frequencies or cell with {wmin,wmax}
%       -LineSpec: String with line style, marker symbol, and color. See <a href="matlab:doc plot">plot</a>
%
% Output Arguments:      
%       -s:     vector of singular values of complex frequency response
%       -omega: vector of complex frequencies
%
% Examples:
%       The following code plots the singular values of the frequency 
%       response of the benchmark 'CDplayer' (SSS, MIMO):
%
%> load CDplayer.mat; sys=sss(A,B,C);
%> figure; sigma(sys);
%
% See Also:
%       bode, freqresp, bodemag, bodeplot
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
%               Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
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

if nargout>0 && length(varargin)>1
    error('The "sigma" command operates on a single model when used with output arguments.');
end

%% Get frd object
for i=1:length(varargin)
    if isa(varargin{i},'sss')
        varargin{i} = frd(varargin{i}, omega);
    end
end

%% Call ss/sigma
if nargout>0 %no plot
    [s, omega]=sigma(varargin{1});
else
    sigma(varargin{:});
end
