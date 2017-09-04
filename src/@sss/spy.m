function varargout = spy(varargin)
% SPY - Plot sparsity pattern of sss system
% 
% Syntax:
%               SPY(sys)
%               SPY(sys,name)
% 
% Description:
%       This function plots the sparsity pattern of the E and A matrices of
%       the sparse state-space system sys into a new figure. 
%
%       It is possible to pass a name via a second optional argument or
%       receive the figure handle as an output. If no name is passed, then
%       the plot title is set to sys.Name
%
% Input Arguments:
%       *Required Input Arguments:*
%		-sys:  sparse state space (sss)-object
%       *Optional Input Arguments:*
%       -name: Plot title
%       
%
% Examples:
%		The following code plots the sparsity pattern of the benchmark
%		'SpiralInductorPeec' (DSSS, SISO):
%
%> sys = sss('SpiralInductorPeec.mat'); 
%> figure; spy(sys,'Peec inductor');
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
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  17 Mar 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

[varargout{1:nargout}] = sssFunc.spy(varargin{:});