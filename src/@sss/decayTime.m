function varargout = decayTime(varargin)
% DECAYTIME - Computes the time period in which a sparse LTI system levels off
%
% Syntax:
%       DECAYTIME(sys)
%       tmax = DECAYTIME(sys)
%
% Description:
%       tmax = DECAYTIME(sys) computes the time tmax in which the sparse LTI
%       system sys levels off. This is done by taking the slowest pole among
%       the dominant ones and then by computing the time in which the
%       slowest pole decays to 1% of its maximum amplitude.
%
%       If sys is not stable (i.e. there exists at least one pole whose
%       real part is >0), then the decay time is set to tmax=NaN and a 
%       warning is displayed.
%
% Input Arguments:
%       -sys: an sss-object containing the LTI system
%
% Output Arguments:
%       -tmax: time after which the system has settled down
%
% Examples:
%       This code computes the decay time of the benchmark 'building':
%
%> load building; 
%> sys = sss(A,B,C);
%> tmax = decayTime(sys)
%
%       You can visualize the meaning of decay time in a step response plot:
%
%> step(sys);
%> hold on; plot(tmax*[1,1],[-8,8]*1e-4,'r');
%
% See Also:
%       residue, step, impulse
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
% Authors:      Heiko Panzer, Sylvia Cremer, Maria Cruz Varona, 
%               Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  19 Jan 2016
% Copyright (c) 2015, 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

[varargout{1:nargout}] = sssFunc.decayTime(varargin{:});

