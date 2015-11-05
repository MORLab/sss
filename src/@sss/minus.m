function diff = minus(sys1, sys2)
% diff - Computes difference of two LTI systems: u-->(sys1-sys2)-->y
% 
% Syntax:
%       diff = minus(sys1, sys2)
%       diff = sys1-sys2 
%
% Description:
%       Computes difference of two LTI systems: u-->(sys1-sys2)-->y
%
% Input Arguments:       
%       -sys1: minuend sss-object
%       -sys2: subtrahend sss-object
%
% Output Arguments:      
%       -sss-object representing sys1-sys2
%
% Examples:
%
% See Also:
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
% Authors:      Heiko Panzer
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if sys1.n == 0
    diff = sss(sys2.A, sys2.B, -sys2.C, sys2.D, sys2.E);
    return
end
if sys2.n == 0
    diff = sss(sys1.A, sys1.B, sys1.C, sys1.D, sys1.E);
    return
end
if sys1.p ~= sys2.p
    error('sys1 and sys2 must have same number of inputs.')
end
if sys1.m ~= sys2.m
    error('sys1 and sys2 must have same number of outputs.')
end

diff = sss([sys1.A sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.A], ...
          [sys1.B; sys2.B], ...
          [sys1.C, -sys2.C], ...
          sys1.D - sys2.D, ...
          [sys1.E sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.E]);
