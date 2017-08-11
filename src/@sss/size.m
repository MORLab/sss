function [varargout] = size(varargin)
% SIZE - Computes the size of a sparse LTI system (sss)
%
% Syntax:
%       size(sys);
%       sizeSys = size(sys);
%       p = size(sys,1);
%       m = size(sys,2);
%       
% 
% Description:
%       size(sys) computes the size of sys (number of outputs, number of inputs
%       and number of states) and displays them on the Command Window.
%
%       sizeSys = size(sys) computes the size of sys (number of outputs and number
%       of inputs) and stores it in the 1x2 vector sizeSys. This means:
%       p=sizeSys(1) and m=sizeSys(2).
%
% Input Arguments:        
%       -sys: sparse state space (sss)-object
%
% Output Arguments:
%       -p:  output dimension
%       -m:  input dimension
%
% Examples:
%       The following code computes the size of the benchmark model
%       'rail_1357' (DSSS, MIMO):
%
%> load rail_1357.mat
%> sys=sss(A,B,C,[],E);
%> sizeRail_1357=size(sys); %Vector with number of outputs and inputs
%> size(sys); %displaying the information in the Command Window
%
%       The 'sss' class has this functionality also implemented. When you
%       want to get the number of output, input and state variables of a
%       system, you just type:
%
%> p = sys.p %Number of outputs
%> m = sys.m %Number of inputs
%> n = sys.n %Number of state variables
%
% See Also:
%       sss, spy, issd
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
% Authors:      Thomas Emmert, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

sys = varargin{1};

if nargin==1
    if nargout==0
        disp(['Sparse state space model with ',num2str(sys.p),' outputs, ',num2str(sys.m),' inputs, and ',num2str(sys.n),' states.'])
        return
    end
    
    varargout{1} = [size(sys.c,1) size(sys.b,2)];

elseif nargin==2
    if varargin{2} == 1
        varargout{1} = size(sys.c,1);
    else
        varargout{1} = size(sys.b,2);
    end
end