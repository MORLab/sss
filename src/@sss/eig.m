function varargout = eig(varargin)
% EIG - Compute eigenvalues and eigenvectors of a sparse state-space model
%
% Syntax:
%       EIG(sys)
%       e = EIG(sys)
%       [V,D] = EIG(sys)
%       [V,D,W] = EIG(sys)
%       [V,D,W] = EIG(sys,varargin)
%
% Description:
%       e = eig(sys) returns a column vector e containing the eigenvalues of 
%       the sparse state-space system sys. If sys is descriptor (sys.E ~= I, 
%       sys.E invertible), then the column vector e contains the generalized
%       eigenvalues of the pencil (sys.A, sys.E).
%
%       [V,D] = eig(sys) returns diagonal matrix D of eigenvalues and matrix V 
%       whose columns are the corresponding right eigenvectors, so that 
%       sys.A*V = sys.E*V*D.
%
%       [V,D,W] = eig(sys) also returns full matrix W whose columns are the
%       corresponding left eigenvectors, so that W'*sys.A = sys.E*D*W'.
%
%       //Note: this function performs dense computations and is hence
%       suitable only for mid-sized problems. Use <a href="matlab:open eigshelp.html">eigs</a> instead if you
%       are interested only in some eigenmodes.
%
% Input Arguments:
%       -sys: an sss-object containing the LTI system
%
% Output Arguments:
%       -e: eigenvalues returned as column vector
%       -V: right eigenvectors (square matrix)
%       -D: diagonal matrix of eigenvalues
%       -W: left eigenvectors (square matrix)
%
% Examples:
%       Load the benchmark 'building' (SSS, SISO) and compute the eigenvalues.
%
%> load building.mat
%> sys = sss(A,B,C);
%> e = eig(sys);
%
%       Load the benchmark 'rail_1357' (DSSS, MIMO) and compute the
%       generalized eigenvalues as well as the right and left eigenvectors
%       of the pair (sys.A, sys.E).
%
%> load rail_1357.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m),E);
%> [V,D,W] = eig(sys);
%
% See Also:
%       matfun/eig, eigs
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
% Authors:      Heiko Panzer, Thomas Emmert (emmert@tfd.mw.tum.de), Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  29 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

[varargout{1:nargout}] = sssFunc.eig(varargin{:});
