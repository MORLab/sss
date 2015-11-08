function [varargout]  = size(varargin)
% SIZE - Computes the size of a sparse LTI system (sss)
%
% Syntax:
%       [p, m] = size(sys);
% 
% Description:
%       Computes the size of a sparse LTI system (sss)
%
% Input Arguments:        
%       -sys: sparse state space (sss)-object
%
% Output Arguments:
%       -p:  output dimension
%       -m:  input dimension
%
% Examples:
%       TODO
%
% See Also:
%       TODO
%
% References:
%       TODO
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Chair of Thermofluid Dynamics, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Thomas Emmert, Maria Cruz Varona
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
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