function varargout = minus(varargin)
% MINUS - Computes difference of two sparse LTI systems.
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
%       -sys1: minuend sss-object
%       -sys2: subtrahend sss-object
%
% Output Arguments:      
%       -diff: sss-object representing sys1-sys2
%
% Examples:
%       In this example the 'building' model will be reduced using the build-in
%       |balancmr| function. Note that for the reduction of large-scale models,
%       we recommend using the <https://www.rt.mw.tum.de/?sssMOR sssMOR> toolbox.
%       we acts directly on |sss| models and exploits sparsity. 
%
%> load building.mat, sys=sss(A,B,C);
%> sysr=sss(balancmr(ss(sys),12)); %reduced order: 12
%
%       The reduced model is compared to the original in a bode magnitude
%       plot. We use the |minus| function to compute the error model |syse|.
%
%> syse=minus(sys,sysr); %syse = sys - sysr
%> figure; bodemag(sys,sysr,'--r',syse,'--g')
%> legend('original','reduced','error')
%
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

[varargout{1:nargout}] = sssFunc.minus(varargin{:});
