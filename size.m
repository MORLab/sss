function [varargout]  = size(varargin)
% Computes the size of a sparse LTI system (sss)
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% [p, m]  = size(sys);
% Input:        * sys: sparse state space (sss)-object
% Output:       * p:  output dimension
%               * m:  input dimension
% ------------------------------------------------------------------
% Authors:      Thomas Emmert (emmert@tfd.mw.tum.de)
%               Maria Cruz Varona (maria.cruz@tum.de)
% Last Change:  22 Sep 2015
% ------------------------------------------------------------------

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