function  varargout = impulse(varargin)
% IMPULSE - Computes and/or plots the impulse response of a sparse LTI system
%
% Syntax:
%   IMPULSE(sys)
%   IMPULSE(sys,t)
%   IMPULSE(sys,Tfinal)
%   IMPULSE(sys1, sys2, ..., t)
%   IMPULSE(sys1, sys2, ..., Tfinal)
%   IMPULSE(sys1,'-r',sys2,'--k',t)
%   IMPULSE(sys1,'-r',sys2,'--k',Tfinal)
%   [h, t] = IMPULSE(sys,t)
%   [h, t] = IMPULSE(sys,Tfinal)
%   [h, t] = IMPULSE(sys,...,Opts)
%   TF     = IMPULSE(sys,...,struct('tf',true))
%   [TF,h,t] = IMPULSE(sys,...,struct('tf',true))
%
% Description:
%       Computes and/or plots the impulse response of a sparse LTI system
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -t:     vector of time values to plot at
%       -Tfinal: end time of impulse response
%       -Opts:  structure with execution parameters
%			-.odeset:  settings of ODE solver
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
%           -.tLin: uniformly spaced time vector
%                       [{0} / 1]
%
% Outputs:
%       -h, t: vectors containing impulse response and time vector
%       -TF: discrete time |tf| object of step response
%
% Examples:
%       The following code computes the impulse response of the benchmark
%       'building' (SSS, SISO) and compares it with the MATLAB built-in function:
%
%> load beam.mat
%> sys=sss(A,B,C);    %sparse state-space (sss)
%> figure; impulse(sys); 
%> hold on; impulse(ss(sys)); %compare to built-in
%> legend('ss/impulse','sss/impulse');
%
% See Also:
%       residue, ss/impulse, step
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

[varargout{1:nargout}] = sss.impulse(varargin{:});
