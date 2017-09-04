function [X,rNormVec] = iterativeRefinement(A,B,X,Opts)
% ITERATIVEREFINEMENT - Use iterative LSE solver to refine LSE solution
% 
% Syntax:
%		X                   = ITERATIVEREFINEMENT(A,B,X)
%       [X,rNormVec]        = ITERATIVEREFINEMENT(A,B,X)
%		[...]               = ITERATIVEREFINEMENT(A,B,X,Opts)
% 
% Description:
%       This function takes as input the approximate solution X to the
%       linear system A*X = B and tries to improve its accuracy. B is also
%       allowed to be a matrix, so that X is the solution of a set of lse
%       with common A matrix.
%
%       This can be done either by direct methods (compare [1,2]) or 
%       by feeding the data to an iterative solver (so far, only cgs is 
%       implemented). In both cases, the algorihm stops when the desired
%       accuracy, measured as the relative norm of the residual
%       r = B-A*X is achieved.
%
%       This can be useful when subspecting that a solution X obtained
%       through direct methods (e.g. lu or chol) may be inaccurate due to
%       roundoff errors.
%
% Input Arguments:
%		*Required Input Arguments:*
%		-A:             Left-hand side matrix of the lse
%		-B:             Right-hand side matrix of the lse
%       -X:             (inaccurate) Solution to the lse
%
%		*Optional Input Arguments:*
%		-Opts:	 		A structure containing following fields
%			-.method:  	refinement method to be executed;
% 						[{wilkinson} / cgs ]
%			-.tol:  	desired tolerace for the residual norm;
% 						[{1e-15}]
%			-.maxiter:  maximum number of iterations;
% 						[{1e2}]
%           -.L/U:      LU factors
% 						[{[ ]}]
%           -.P/Q/S:    Permutation matrices from the sparse LU
% 						[{speye(size(A))}]
%
% Output Arguments:
%       -X:             The solution to the lse satisfying the desired accuracy 
%       -rNormVec:      Vector of relative norm of R through the iterations
%                       //Note: if Opts.method='cgs' and size(B,2)>1, then only the rNormVec for the last column of X is returned
%
% Examples:
%       In this example, the solution of a linear system of equations is
%       achieved through \ and is subsequently perturbed to decrease
%       accuracy. Through a call to ITERATIVEREFINEMENT, the accurate
%       solution is regained
%
%>sys = sss('CDplayer');
%>A = sys.A; B = sys.B; X = A\B; 
%>R = A*X-B; rNormOrig = norm(R,'fro')/norm(B,'fro')
%>
%>X = X + 1e6*eps*X; R = A*X-B; rNormPert = norm(R,'fro')/norm(B,'fro')
%>
%>X = iterativeRefinement(A,B,X); 
%>R = A*X-B; rNormRef = norm(R,'fro')/norm(B,'fro')
%
%
% See Also: 
%		lu, chol, mldivide, cgs
%
% References:
%		* *[1] Wilkinson (1963)*, Rounding Errors in Algebraic Processes
%		* *[2] Moler (1967)*, Iterative Refinement in Floating Point
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  06 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Parsing
Def.method  = 'wilkinson';
Def.tol     = 1e-15;
Def.maxiter = 1e2;
% from: [L,U,P,Q,S] = lu(A) where P*(S\A)*Q = L*U, hence X=Q*(U\(L\(P*(S\B))))
Def.L   = []; 
Def.U   = [];
Def.P   = speye(size(A));
Def.Q   = speye(size(A));
Def.S   = speye(size(A));

if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end    

%% Refinement

switch Opts.method
    case 'wilkinson'
        bNorm = norm(B,'fro');
        k = 0; rNormVec = zeros(1,Opts.maxiter);
        
        R = B-A*X; %residual
        rNorm = norm(R,'fro')/bNorm;
        
        %   Compute LU if not available, provided we need to improve R
        if rNorm > Opts.tol && (isempty(Opts.L) || ~isempty(Opts.U))
            [Opts.L, Opts.U, Opts.P, Opts.Q, Opts.S] = lu(A);
        end
        
        %   Refine!
        while rNorm > Opts.tol && k <= Opts.maxiter
            k = k+1;

            D = Opts.Q*(Opts.U\(Opts.L\(Opts.P*(Opts.S\R))));
            X = X + D;      
            
            R = B-A*X;
            rNorm = norm(R,'fro')/bNorm;
            rNormVec(k) = rNorm;
        end
        
        if k <= Opts.maxiter, rNormVec(k+1:end) = []; end
            
    case 'cgs'
        warning('off','MATLAB:cgs:tooSmallTolerance');
        for iCol = 1:size(X,2)
            %Solve LSE with desired accuracy
            [temp,exitFlag,rNormVec] = cgs(A,B(:,iCol),Opts.tol,Opts.maxiter,Opts.L,Opts.U,X(:,iCol));
            %replace column if refinement worked
            if exitFlag == 0,
                X(:,iCol) = temp;
            else
                warning('sss:iterativeRefinement:cgsNoConvergence','iterative refinement did not work.');
            end
        end
        warning('on','MATLAB:cgs:tooSmallTolerance');
end



