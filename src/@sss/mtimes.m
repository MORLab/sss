function prod = mtimes(sys1, sys2)
% MTIMES
%
% Syntax:
%       prod = MTIMES(sys1, sys2)
%       prod = sys1*sys2
%
% Description:
%       mtimes computes the product of two LTI systems: u --> sys2 --> sys1
%       --> y, i.e., the output of sys2 is connected directly to the input
%       of sys1.
%
% Input Arguments:
%       sys1, sys2: sss-object containing the LTI system
% 
% Outputs Arguments:
%       prod: sss-object representing sys1*sys2
%
% See Also: 
%
% Examples:
%
% References:
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%------------------------------------------------------------------
% Authors:      Heiko Panzer
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if (sys1.n==0 && isempty(sys1.D)) || (sys2.n==0 && isempty(sys2.D))
    prod = sss([], [], []);
    return
end
if sys2.p ~= sys1.m
    error('Number of outputs of sys1 does not match number of inputs of sys2.')
end

prod = sss([sys1.A sys1.B*sys2.C; sparse(sys2.n,sys1.n) sys2.A], ...
          [sys1.B*sys2.D; sys2.B], ...
          [sys1.C, sys1.D*sys2.C], ...
          sys1.D*sys2.D, ...
          [sys1.E sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.E]);
