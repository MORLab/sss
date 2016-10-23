function  varargout = step(varargin)
% STEP - Computes and/or plots the step response of a sparse LTI system
%
% Syntax:
%   STEP(sys)
%   STEP(sys,t)
%   STEP(sys,Tfinal)
%   STEP(sys1, sys2, ..., t)
%   STEP(sys1, sys2, ..., Tfinal)
%   STEP(sys1,'-r',sys2,'--k',t);
%   STEP(sys1,'-r',sys2,'--k',Tfinal)
%   [h, t] = STEP(sys)
%   [h, t] = STEP(sys, t)
%   [h, t] = STEP(sys, Tfinal)
%   [h, t] = STEP(sys, ..., Opts)
%   TF = STEP(sys,...,struct('tf',true))
%   [TF,h,t] = STEP(sys,...,struct('tf',true))
%
% Description:
%       step(sys) plots the step response of the sparse LTI system sys
%
%       [h, t] = step(sys, t) computes the step response of the sparse LTI
%       system sys and returns the vectors h and t with the response and
%       the time series, respectively.
%
%       TF = STEP(sys,struct('tf',true)) returns a discrete time |tf| object 
%       of the FIR-filter with same discrete step response as sys.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -t:     vector of time values to plot at
%       -Tfinal: end time of step response
%       -Opts:  structure with execution parameters
%			-.odeset:  odeset Settings of ODE solver
%           -.tolOutput: Terminate if norm(y_-yFinal)/norm(yFinal)<tolOutput with yFinal = C*xFinal+D;
%						[1e-3 / positive float]
%           -.tolState: Terminate if norm(x-xFinal)/norm(xFinal)<tolState with xFinal = -(A\B);
%						[1e-3 / positive float]
%           -.tf: return tf object
%                       [{0} / 1]
%           -.ode: ode solver;
%                       [{'ode45'} / 'ode113' / 'ode15s' / 'ode23'] 
%           -.tsMin: minimum sample time if no time vector is specified
%                       [{0} / positive float]
%           -.htCell: return ode output as cell with irregularly spaced t
%                       [{0} / 1]
%           -.tLin: uniformly spaced time vector
%                       [{0} / 1]
% 
% Output Arguments:
%       -h, t: vectors containing step response and time vector
%       -TF: discrete time tf object of step response
%
% Examples:
%       The following code plots the step response of the benchmark
%       'building' (SSS, SISO):
%
%> load building.mat; sys=sss(A,B,C);
%> step(sys);
%
% See Also:
%       residue, impulse
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
% Authors:      Heiko Panzer, Sylvia Cremer, Jorge Luiz Moreira Silva
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  14 Jun 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

[varargout{1:nargout}] = sss.step(varargin{:});
