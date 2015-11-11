function ispd = ispd(A)
% ISPD - Determines if a matrix is positive definite 
%
% Syntax:
%       ISPD(A)
%       ispd = ISPD(A)
%
% Description:
%       Positive (or negative) definiteness is an important property of
%       matrices that can be exploited if given. Therefore, it is often
%       desirable to determine whether a matrix has this property or not.
%
%       This function tries a Cholesky factorization of A. If it succedes, 
%       then the matrix is positive definite.
%
% Input Arguments:
%       -A: matrix which have to be checked for positive definitness
%
% Output Arguments:
%       -ispd: a boolean value (1=true, 0=false)
%
% Examples:
%       TODO
%
% See Also:
%       chol, issd
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  11 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Computation
[~,p] = chol(A); 
ispd = ~p>0;

%%  Print result if no ouput was defined
if nargout == 0
    if p>0
        fprintf('Matrix is NOT positive definite \n');
    else
        fprintf('Matrix IS positive definite \n');
    end
end
