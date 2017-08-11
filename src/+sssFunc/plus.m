function sum = plus(sys1, sys2)
% PLUS - Computes sum of two sparse LTI systems.
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
%       -sys1, sys2: summand sss-object
%
% Output Arguments:      
%       -sum: sss-object representing sys1+sys2
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
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  08 Apr 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
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


% define system size, because sys.m, sys.n and sys.p are not defined for
% ss-objects
sys1n = size(sys1.A,1);
sys2n = size(sys2.A,1);
sys1m = size(sys1.B,2);
sys2m = size(sys2.B,2);
sys1p = size(sys1.C,1);
sys2p = size(sys2.C,1);

% change sys.E = [] to sys.E = eye(n)
if isempty(sys1.E) sys1.E=sparse(eye(sys1n)); end
if isempty(sys2.E) 
    if isa(sys2,'sss')
        sys2.E=sparse(eye(sys2n));
    else
        sys2.E=eye(sys2n);
    end
end

if sys1n == 0, sum = sys2; return, end
if sys2n == 0, sum = sys1; return, end
if sys1m ~= sys2m
    error('sys1 and sys2 must have same number of inputs.')
end
if sys1p ~= sys2p
    error('sys1 and sys2 must have same number of outputs.')
end

% Create the sum object
if isa(sys1,'sss') || isa(sys2,'sss')
    sum = sss([sys1.A sparse(sys1n,sys2n); sparse(sys2n,sys1n) sys2.A], ...
          [sys1.B; sys2.B], ...
          [sys1.C, sys2.C], ...
          sys1.D + sys2.D, ...
          [sys1.E sparse(sys1n,sys2n); sparse(sys2n,sys1n) sys2.E]);
else %ssRed
    sum = ssRed([sys1.A zeros(sys1n,sys2n); zeros(sys2n,sys1n) sys2.A], ...
          [sys1.B; sys2.B], ...
          [sys1.C, sys2.C], ...
          sys1.D + sys2.D, ...
          [sys1.E zeros(sys1n,sys2n); zeros(sys2n,sys1n) sys2.E]);
end
