function varargout = plus(varargin)
% PLUS - Computes sum of two sparse LTI systems.
% 
% Syntax:
%       sum = PLUS(sys1, sys2)
%       sum = sys1+sys2
% 
% Description:
%       PLUS gives as output the system sum, which is the combination of two 
%       different systems sys1 and sys2 that have the same number of inputs/outputs. 
%       The output of the system sum will be equivalent to the sum of the 
%       outputs from sys1 and sys2, considering that they are excited with 
%       the same input u (u-->(sys1+sys2)-->y).
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
% Authors:      Heiko Panzer
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

[varargout{1:nargout}] = sssFunc.plus(varargin{:});
