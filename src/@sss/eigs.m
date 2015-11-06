function varargout = eigs(sys, varargin)
% EIGS - compute eigenvalues of the sparse state space system using sparse matrices.
%
% Syntax:
%        D = eigs(sys)
%       [V, D, flag] = eigs(sys)
%       eigs(sys,k)
%       eigs(sys,k,sigma)
%       eigs(sys,k,sigma,opts)
%
% Description:
%       Compute eigenvalues of the sparse state space system using sparse
%       matrices.
% 
% Input Arguments:
%       *required Input Arguments:*
%       -sys:   An sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -k,sigma,opts:  Output and computation options,
%        see <a href="matlab:help eigs">eigs</a>
%
% Output:       
%       -V:     Right eigenvectors
%       -D:     Diagonal matrix of eigenvalues
%       -flag:  Convergence flag
%
% Examples:
%
% See Also: 
%        eigs
%
% References:
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Chair of Thermofluid Dynamics, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Thomas Emmert
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

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