function isspd = isspd(A)
% ISSPD - Determines if a matrix is symmetric positive definite 
%
% Syntax:
%       ISSPD(A)
%       isspd = ISSPD(A)
%
% Description:
%       Symmetric positive (or negative) definiteness is an important property of
%       matrices that can be exploited if given. Therefore, it is often
%       desirable to determine whether a matrix has this property or not.
%
%       This function first checks if the matrix A is symmetric. If the 
%       matrix is nearly symmetric, then the function symmetrizes it by 
%       computing the symmetric part via (A + A.')/2, where norm(A-Asymm')
%       Then, it tries a Cholesky factorization of A. If it succeeds, then 
%       the matrix is symmetric positive definite.
%
% Input Arguments:
%       *Required Input Arguments:*
%           -A: matrix which have to be checked for positive definitness
%
% Output Arguments:
%       -isspd: a boolean value (1=true, 0=false)
%
%
% See Also:
%       chol, ispd, issd
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
% Authors:      Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  26 Oct 2018
% Copyright (c) 2015-2018 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Computation
if issymmetric(A) && norm(A-A.','fro') < 1e-8 % ||, &&
    fprintf('Matrix is symmetric \n');
    [~,p] = chol(A); 
    isspd = (~p)>0;
else
    Asymm = (A + A.')/2; % compute symmetric part of A
    if issymmetric(Asymm) && norm(full(A-Asymm)) < 1e-8
        fprintf('Matrix is nearly symmetric and was symmetrized by computing the symmetric part (A + A.'')/2 \n');
        [~,p] = chol(Asymm);
        isspd = (~p)>0;
    else
        error('Matrix is not symmetric \n');
    end
end

%%  Print result if no ouput was defined
if nargout == 0
    if p>0
        fprintf('Matrix is NOT positive definite \n');
    else
        fprintf('Matrix IS positive definite \n');
    end
end
