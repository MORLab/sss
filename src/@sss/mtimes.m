function sys1 = mtimes(sys1, sys2)
% MTIMES - Computes the product of two LTI systems: u-->sys2-->sys1-->y
%
% Syntax:
%       prod = MTIMES(sys1, sys2)
%       prod = sys1*sys2
%
% Description:
%       mtimes computes the product of two LTI systems: prod = sys1*sys2,
%       i.e., the output of sys2 is connected directly to the input
%       of sys1.
%
% Input Arguments:
%       sys1, sys2: sss-object containing the LTI system
% 
% Outputs Arguments:
%       prod: sss-object representing sys1*sys2
%
% Examples:
%> load building.mat
%> sysBuilding=sss(A,B,C);
%> size(sysBuilding)
%> sysRandom=rss(10); sysRandomSparse=sss(sysRandom);
%> size(sysRandomSparse);
%> prodSystems=mtimes(sysBuilding,sysRandomSparse);
%> size(prodSystems);
%> prod = sysBuilding*sysRandomSparse; %for comparison reasons
%> prod.A == prodSystems.A; %for comparison reasons
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
%                   -> <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Thomas Emmert
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


%% Convert numeric gain to sss
if isnumeric(sys1)
    D = sparse(sys1);
    if isvector(D)
        % Apply static gain to all channels
        D = speye(sys2.p)*diag(D);
    end
    % Give all properties to the returned model
    sys1 = sys2.clear;
    sys1.D = D;
    sys1.y = sys2.y;
end
if isnumeric(sys2)
    D = sparse(sys2);
    if isvector(D)
        % Apply static gain to all channels
        D = speye(sys1.m)*diag(D);
    end
    sys2 = sss(D);
    sys2.u = sys1.u;
end

if sys1.m ~= sys2.p
    error('Number of inputs of sys1 does not match number of outputs of sys2.')
end

% Store input and output names of the resulting system
u = sys2.u;

A = [sys1.A sys1.B*sys2.C; sparse(sys2.n,sys1.n) sys2.A];
B = [sys1.B*sys2.D; sys2.B];
C = [sys1.C, sys1.D*sys2.C];
D = sys1.D*sys2.D;
E = [sys1.E sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.E];

sys1.A=A; sys1.B=B; sys1.C=C; sys1.D=D; sys1.E=E;

sys1.u = u;

end