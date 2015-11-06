function varargout = eig(sys, varargin)
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
%       eigenvalues of the pencil (sys.A,sys.E).
%
%       [V,D] = eig(sys) returns diagonal matrix D of eigenvalues and matrix V 
%       whose columns are the corresponding right eigenvectors, so that 
%       sys.A*V = sys.E*V*D.
%
%       [V,D,W] = eig(sys) also returns full matrix W whose columns are the
%       corresponding left eigenvectors, so that W'*sys.A = sys.E*D*W'.
%
%       If the dimension of sys is too large (for instance, sys.n>10000),
%       then it is preferable to use eigs(sys) instead.
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
%       Load the benchmark "build" (SSS,SISO) and compute the eigenvalues.
%
%> load build.mat
%> sys = sss(A,B,C)
%> e = eig(sys);
%
%       Load the benchmark "PEEC_MTLn1600" (DSS,MIMO) and compute the
%       generalized eigenvalues as well as the right and left eigenvectors
%       of the pair (sys.A,sys.E).
%
%> load PEEC_MTLn1600.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m),E)
%> [V,D,W] = eig(sys);
%
% See Also:
%       ss/eig, eigs
%
% ------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
% ------------------------------------------------------------------
% Authors:      Heiko Panzer, Thomas Emmert (emmert@tfd.mw.tum.de),
%               Maria Cruz Varona
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  29 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

if sys.isBig
    warning(['System order is very large: ',num2str(sys.n),'. You may want to try eigs(sys) instead.'])
end

if sys.isDescriptor
    if nargout==1||nargout==0
        [varargout{1}] = eig(full(sys.a), full(sys.e),varargin{:});
    elseif nargout == 2
        [varargout{1}, varargout{2}] = eig(full(sys.a), full(sys.e),varargin{:});
    elseif nargout == 3
        [varargout{1}, varargout{2}, varargout{3}]  = eig(full(sys.a), full(sys.e),varargin{:});
    end
else
    if nargout==1||nargout==0
        [varargout{1}] = eig(full(sys.a),varargin{:});
    elseif nargout == 2
        [varargout{1}, varargout{2}] = eig(full(sys.a),varargin{:});
    elseif nargout == 3
        [varargout{1}, varargout{2}, varargout{3}]  = eig(full(sys.a),varargin{:});
    end
end

for i = 1:nargout
    varargout{i} = varargout{i};
end
% Store poles for future computations
if nargout>1
    sys.poles = varargout{2};
else
    sys.poles = varargout{1};
end

end