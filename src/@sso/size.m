function [varargout] = size(varargin)
% SIZE - Computes the size of a sparse second-order system (sso)
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
%       -sys: sparse second-order (sso)-object
%
% Output Arguments:
%       -p:  output dimension
%       -m:  input dimension
%
% Examples:
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

sys = varargin{1};

if nargin==1
    if nargout==0
        disp(['Sparse model model with ',num2str(sys.p),' output(s), ',num2str(sys.m),' input(s), and ',num2str(sys.n),' states.'])
        return
    end
    
    varargout{1} = [size(sys.Cp,1) size(sys.B,2)];

elseif nargin==2
    if varargin{2} == 1
        varargout{1} = size(sys.Cp,1);
    else
        varargout{1} = size(sys.B,2);
    end
end