function varargout = issd(varargin)
% ISSD - Check strict dissipativity of sparse LTI system
%
% Syntax:
%       ISSD(sys)
%       issd = ISSD(sys)
%       [issd,numericalAbscissa] = ISSD(sys)
%
% Description:
%       This function determines whether the LTI, sss system |sys| is 
%       given in a strictly dissipative realization, satisfying
%       $$ E = E^T > 0, \quad A+A^T<0. $$  
%       The computations are meant 
%       to avoid dense operations. However, whenever this is not possible, 
%       a warning is issued.
%
%       If no output is defined, then the result is printed on the screen.
%       Depending on the number of ouputs defined, the function can return.
%
%       NaN is returned either when the computation was not possible, or 
%       when the numerical abscissa is zero. In the latter case, the system
%       might be stable (in the sense of Lyapunov) or unstable if the
%       multiplicity of the eigenvalues at the origin is greater than one.
%
% Input Arguments:
%       -sys: sss-object containing the LTI system
%
% Output Arguments:
%       -issd: a boolean value (1=true, 0=false, NaN= dissipative but not
%              strictly).
%       -numericalAbscissa: the largest eigenvalue of A+A'.
%
% Examples:
%       The following code checks the strictly dissipativity of the
%       benchmark 'CDplayer' (SSS, MIMO):
%> clear all;
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> [issd, numericalAbscissa]=issd(sys)
%	
% See Also:
%       ispd, eigs, chol, sparse
%
% References:
%       * *[1] Panzer (2014)*, Model order reduction by Krylov subspace methods
%       with Global Error Bounds and Automatic Choice of Parameters
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
% Authors:      Heiko Panzer, Alessandro Castagnotto, Maria Cruz Varona 
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  31 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

[varargout{1:nargout}] = sssFunc.issd(varargin{:});

