function [y,x_,k,index] = simUpdate(y,x_,k,index,x,u,i,m,C,D)  
% SIMUPDATE - Update for simulation functions
%
% Syntax:
%       y              = SIMUPDATE(y,x_,k,index,x,u,i,m,C,D)
%       [y,x_,k]       = SIMUPDATE(y,x_,k,index,x,u,i,m,C,D)
%       [y,x_,k,index] = SIMUPDATE(y,x_,k,index,x,u,i,m,C,D)
%
% Description:
%       Auxiliary function for the update of state and ouput after each 
%       computation step of the simulation functions.
%
% Input Arguments:
%       -y:         output vector
%       -x_:        matrix of state vectors
%       -k,index,i: loop variables
%       -u, x:      input matrix and and current state vector
%       -m:         see simulation functions (e.g. simRK4 for more info)
%       -C,D:       matrices of state-space model
%
% Output Arguments:
%       -y:         output vector
%       -x_:        matrix of state vectors
%       -k,index:   loop variables
%
% See Also:
%       sim, simForwardEuler, simRK4, simInit
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
% Last Change:  05 Sep 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% Update output and state
y(:,i) = C*x + D*u(i,:)';
if ~isempty(x_)
    if mod(i,m) == 0
        k       = k+1;            
        x_(:,k) = x;
        index   = [index i];
    end
end