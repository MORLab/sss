function varargout = plus(varargin)
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


[varargout{1:nargout}] = sssFunc.plus(varargin{:});
