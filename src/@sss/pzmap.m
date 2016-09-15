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
        if varargin{i}.isDae
            error('pzmap does not work with DAE systems yet.');
        end
        varargin{i} = zpk(varargin{i}, k, struct('zpk',true));
    end
end

if nargout
    [varargout{1}, varargout{2}]=pzmap(varargin{:});
else
    pzmap(varargin{:});
end
end


