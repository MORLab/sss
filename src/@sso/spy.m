function varargout = spy(varargin)
% SPY - Plot sparsity pattern of sso system
% 
% Syntax:
%               SPY(sys)
%               SPY(sys,name)
% 
% Description:
%       This function plots the sparsity pattern of the M, D and K matrices of
%       the sparse sectond-order system sys into a new figure. 
%
%       It is possible to pass a name via a second optional argument or
%       receive the figure handle as an output. If no name is passed, then
%       the plot title is set to sys.Name
%
% Input Arguments:
%       *Required Input Arguments:*
%		-sys:  sparse second-order (sso)-object
%       *Optional Input Arguments:*
%       -name: Plot title
%       
%
% Examples:
%
%
% See Also: 
%		spy
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

[varargout{1:nargout}] = ssoFunc.spy(varargin{:});