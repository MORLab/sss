function varargout = lyapchol(varargin)
% LYAPCHOL - Solve Lyapunov equations
%
% Syntax:
%       S				= LYAPCHOL(sys)
%       [S,R]			= LYAPCHOL(sys)
%       [S,R]  		    = LYAPCHOL(sys,Opts)
%
% Description:
%       This function returns the Cholesky factorization X=S*S' of the 
%       solution of the Lyapunov equation A*X+X*A'+B*B'=0 or the generalized 
%       Lyapunov equation A*X*E'+E*X*A'+B*B'=0.
%
%       If the number of output arguments is 2, then the low rank factor 
%       Y = R*R' of the dual (generalized) lyapunov equation 
%       A'*Y*E+E'*Y*A+C'*C=0 is solved as well.
%
%       If the option 'type' is set to 'adi',then a low rank approximation 
%       of the Cholesky (like) factor is performed [1]. If this option is not 
%       specified, then ADI is applied to systems with n>500. The options 
%       'lse', 'rctol' and 'q' only apply to ADI.
%
%       //Note: the definition of the Cholesky factors X = S*S' is
%       different from built-in lyapchol, where X = S'*S. However, our
%       definition is consistent both with standard literature (cp [3]) and
%       the low-rank approximation if R has fewer columns than rows.
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%		*Optional Input Arguments:*
%       -Opts:              a structure containing following options
%           -.method:       select solver for lyapunov equation 
%                           [{'auto'} / 'adi' / 'hammarling' / 'crksm']
%           -.lse:          solve linear system of equations (only ADI)
%                           [{'gauss'} / 'luChol']
%           -.rctol:        tolerance for difference between ADI iterates
%                           [{'1e-12'} / positive float]
%           -.q:            size of Cholesky factor (only ADI)
%                           [{'0'} / positive integer]
%           -.forceOrder:   return order q
%                           [{'false'} / 'true']
%           -.maxiter:      maximum number of iterations (only ADI)
%                           [{300} / positive integer]
%
% Output Arguments:
%       -S:     Cholesky factor X=S*S' of Lyapunov equation A*X*E'+E*X*A'+B*B'=0
%       -R:     Cholesky factor Y=R*R' of Lyapunov equation A'*Y*E+E'*Y*A+C'*C=0
%
% Examples:
%       Compute the Cholesky factors for both Lyapunov equations
%
%> sys = sss('building');
%> [S,R] = lyapchol(sys);
%
%       To compute a single Cholesky factor, use
%
%> S = lyapchol(sys);
%> R = lyapchol(sys');
%
% See Also:
%       solveLse, tbr, norm, numerics/lyapchol
%
% References:
%       * *[1] Penzl (2000)*, LYAPACK - A MATLAB Toolbox for Large Lyapunov
%       and Riccati Equations, Model Reduction Problems, and Linear-Quadratic 
%       Optimal Control Problems.
%       * *[2] Saak, Köhler, Benner (2016)*, M-M.E.S.S. - The Matrix 
%       Equation Sparse Solver Library.
%       * *[3] Golub, Van Loan (1996)*, Matrix computations
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
% Authors:      Alessandro Castagnotto, Lisa Jeschek
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  29 Mar 2017
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

[varargout{1:nargout}] = sssFunc.lyapchol(varargin{:});
