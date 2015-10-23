function varargout = eig(sys, varargin)
% compute all eigenvalues of the sparse state space system
% ------------------------------------------------------------------
% [V, D, W] = eig(sys, varargin)
% Input:        * sys: an sss-object containing the LTI system
% Output:       * V right eigenvectors
%               * D diagonal matrix of eigenvalues
%               * W left eigenvectors
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), 
%               Thomas Emmert (emmert@tfd.mw.tum.de)
% Last Change:  25 Feb 2015
% ------------------------------------------------------------------
% see also: eig

if sys.n>10000
    warning(['System order is very large: ',num2str(sys.n),'. Compute time will be very long'])
end

if any(any(sys.E-speye(size(sys.E))))
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

for i =1:nargout
    varargout{i} = sparse(varargout{i});
end
% Store poles for future computations
if nargout>1
    sys.poles = varargout{2};
else
    sys.poles = varargout{1};
end

end