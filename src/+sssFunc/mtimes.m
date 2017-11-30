function prod = mtimes(sys1, sys2)
% MTIMES - Computes the product of two LTI systems.
%
% Syntax:
%       prod = MTIMES(sys1, sys2)
%       prod = sys1*sys2
%
% Description:
%       mtimes computes the product of two LTI systems: prod = sys1*sys2,
%       i.e., the output of sys2 is connected directly to the input
%       of sys1 accortind to u-->sys2-->sys1-->y.
%
% *Required Input Arguments:*
%		- sys1, sys2:       sss-objects containing the LTI systems
% 
% *Outputs Arguments*:
%       - prod:     sss-object representing sys1*sys2
%
% See Also:
%       plus, minus
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
% Authors:      Heiko Panzer, Thomas Emmert
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


%% Convert numeric gain to sss
if isnumeric(sys1)
    k = sys1;

    prod    = sys2;
    prod.D  = prod.D*k;
    prod.C  = prod.C*k;
    
    return
end
if isnumeric(sys2)
  
    k = sys2;
    
    prod    = sys1;
    prod.D  = prod.D*k;
    prod.B  = prod.B*k;
        
    return
end

% Define system size, because sys.m, sys.n and sys.p are not defined for
% ss-objects
sys1n = size(sys1.A,1);
sys2n = size(sys2.A,1);
sys2p = size(sys2.B,2);
sys1m = size(sys1.C,1);

% change sys.E = [] to sys.E = eye(n)
if isempty(sys1.E) sys1.E=sparse(eye(sys1n)); end
if isempty(sys2.E) 
    if isa(sys2,'sss')
        sys2.E=sparse(eye(sys2n));
    else
        sys2.E=eye(sys2n);
    end
end

if sys1m ~= sys2p
    error('Number of inputs of sys1 does not match number of outputs of sys2.')
end

% Store input and output names of the resulting system
u = sys2.u;

if isa(sys1,'sss') || isa(sys2,'sss')
    prod = sss([sys1.A sys1.B*sys2.C; sparse(sys2n,sys1n) sys2.A], ...
          [sys1.B*sys2.D; sys2.B], ...
          [sys1.C, sys1.D*sys2.C], ...
          sys1.D*sys2.D, ...
          [sys1.E sparse(sys1n,sys2n); sparse(sys2n,sys1n) sys2.E]);
else %ssRed
        prod = ssRed([sys1.A sys1.B*sys2.C; zeros(sys2n,sys1n) sys2.A], ...
          [sys1.B*sys2.D; sys2.B], ...
          [sys1.C, sys1.D*sys2.C], ...
          sys1.D*sys2.D, ...
          [sys1.E zeros(sys1n,sys2n); zeros(sys2n,sys1n) sys2.E]);
end

prod.u = u;

end