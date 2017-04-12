function diff = minus(sys1, sys2)
% MINUS - Computes difference of two sso models.
% 
% Syntax:
%       diff = MINUS(sys1, sys2)
%       diff = sys1-sys2 
%
% Description:
%       diff = MINUS(sys1, sys2) computes the difference of the two LTI
%       systems: diff = sys1-sys2
%
% Input Arguments:       
%       -sys1: minuend sso-object
%       -sys2: subtrahend sso-object
%
% Output Arguments:      
%       -diff: sso-object representing sys1-sys2
%
% See Also:
%       plus, mtimes
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Professur fuer Thermofluiddynamik, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  10 Apr 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if sys1.n == 0
    diff = sys2;
    return
end
if sys2.n == 0
    diff = sys2; diff.B = -diff.B;
    return
end
if sys1.m ~= sys2.m
    error('sys1 and sys2 must have same number of inputs.')
end
if sys1.p ~= sys2.p
    error('sys1 and sys2 must have same number of outputs.')
end

diff = sso([sys1.M sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.M], ...
           [sys1.D sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.D], ...
           [sys1.K sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.K], ...
           [sys1.B;     sys2.B], ...
           [sys1.Cp,   -sys2.Cp], ...
           [sys1.Cv,   -sys2.Cv], ...
           sys1.Df     - sys2.Df);
