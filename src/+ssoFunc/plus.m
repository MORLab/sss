function sum = plus(sys1, sys2)
% PLUS - Computes sum of two sso models.
% 
% Syntax:
%       sum = PLUS(sys1, sys2)
%       sum = sys1+sys2
%       sum = sys1+D
% 
% Description:
%       PLUS gives as output the system sum, which is the combination of two 
%       different systems sys1 and sys2 that have the same number of inputs/outputs. 
%       The output of the system sum will be equivalent to the sum of the 
%       outputs from sys1 and sys2, considering that they are excited with 
%       the same input u (u-->(sys1+sys2)-->y).
%
%       The second argument can also be a feedthrough matrix D with same
%       count of input/outputs as sys1.
%
% Input Arguments:
%       -sys1, sys2: summand sso-object
%
% Output Arguments:      
%       -sum: sso-object representing sys1+sys2
%
%
% See Also:
%       minus, mtimes
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
% Last Change:  11 Apr 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% parse inputs
if isa(sys2,'numeric') %feedthrough matrix
    if sys1.m ~= size(sys2,2)
        error('sys1 and sys2 must have same number of inputs.')
    end
    if sys1.p ~= size(sys2,1)
        error('sys1 and sys2 must have same number of outputs.')
    end
    sum = sys1; sum.D = sys1.D + sys2;
    return
end
if sys1.n == 0
    sum = sys2;
    return
end
if sys2.n == 0
    sum = sys2; sum.B = -sum.B;
    return
end
if sys1.m ~= sys2.m
    error('sys1 and sys2 must have same number of inputs.')
end
if sys1.p ~= sys2.p
    error('sys1 and sys2 must have same number of outputs.')
end

sum = sso([sys1.M sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.M], ...
           [sys1.D sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.D], ...
           [sys1.K sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.K], ...
           [sys1.B;     sys2.B], ...
           [sys1.Cp,   +sys2.Cp], ...
           [sys1.Cv,   +sys2.Cv], ...
           sys1.Df     +sys2.Df);