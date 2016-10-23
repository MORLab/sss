function varargout = lyapchol(varargin)
% LYAPCHOL - Solve Lyapunov equations
%
% Syntax:
%       R				= LYAPCHOL(sys)
%       [R,L]			= LYAPCHOL(sys)
%       [R,L]  		    = LYAPCHOL(sys,Opts)
%
% Description:
%       This function returns the Cholesky factorization X=R'*R of the 
%       solution of the Lyapunov equation A*X+X*A'+B*B'=0.
%
%       If the option 'type' is set to 'adi',then a low rank approximation 
%       of the Cholesky factor [1] is performed. If this option is not 
%       specified, then ADI is applied to systems with n>500. The options 
%       'lse', 'rctol' and 'q' only apply to ADI.
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%		*Optional Input Arguments:*
%       -Opts:              a structure containing following options
%           -.type:         select amongst different tbr algorithms
%                           [{''} / 'adi' / 'builtIn' ]
%           -.lse:          solve linear system of equations (only ADI)
%                           [{'gauss'} / 'luChol']
%           -.rctol:        tolerance for difference between ADI iterates
%                           [{'1e-12'} / positive float]
%           -.q:            size of Cholesky factor (only ADI)
%                           [{'0'} / positive integer]
%
% Output Arguments:
%       -R:     Cholesky factor X=R'*R of lyapunov equation A*X+X*A'+B*B'=0
%       -L:     Cholesky factor X=L'*L of lyapunov equation A'*X+X*A+C'*C=0
%
% Examples:
%       Compute the Cholesky factors for both Lyapunov equations
%
%> sys = loadSss('building');
%> [R,L] = lyapchol(sys);
%
%       To compute a single Cholesky factor, use
%
%> R = lyapchol(sys);
%> L = lyapchol(sys');
%
% See Also:
%       solveLse, tbr, norm, isrk
%
% References:
%       * *[1] Penzl (2000)*, LYAPACK - A MATLAB Toolbox for Large Lyapunov
%       and Riccati Equations, Model Reduction Problems, and Linear-Quadratic 
%       Optimal Control Problems.
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
% Authors:      Alessandro Castagnotto, 
%               Lisa Jeschek
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Aug 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

[varargout{1:nargout}] = sss.lyapchol(varargin{:});
