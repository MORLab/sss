function varargout = minus(varargin)
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

[varargout{1:nargout}] = ssoFunc.minus(varargin{:});
