function varargout = eigs(sys, varargin)
% compute eigenvalues of the sparse state space system using sparse
% matrices.
% ------------------------------------------------------------------
% [V, D, flag] = eig(sys, varargin)
% Input:        * sys: an sss-object containing the LTI system
% Output:       * V right eigenvectors
%               * D diagonal matrix of eigenvalues
%               * flag convergence flag
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Thomas Emmert (emmert@tfd.mw.tum.de)
% Last Change:  25 Feb 2015
% ------------------------------------------------------------------
% see also: eigs

if any(any(sys.E-speye(size(sys.E))))
    if nargout==1||nargout==0
        [varargout{1}] = eigs(sys.a, sys.e, varargin{:});
    else
        [varargout{1}, varargout{2}, varargout{3}]  = eigs(sys.a, sys.e, varargin{:});
    end
else
    if nargout==1||nargout==0
        [varargout{1}] = eigs(sys.a, varargin{:});
    else
        [varargout{1}, varargout{2}, varargout{3}]  = eigs(sys.a, varargin{:});
    end
end

% Store poles for future computations
if nargout>1
    sys.poles = varargout{2};
else
    sys.poles = varargout{1};
end

end