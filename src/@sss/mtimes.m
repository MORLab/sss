function varargout = mtimes(varargin)
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

[varargout{1:nargout}] = sssFunc.mtimes(varargin{:});