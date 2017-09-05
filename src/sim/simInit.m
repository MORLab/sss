function [y,x_,m,k,index,L,U,p] = simInit(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor,oneOutput)
% SIMINIT - Initialization for simulation functions
%
% Syntax:
%       y                      = SIMINIT(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor,oneOutput)
%       [y,x_,m,k]             = SIMINIT(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor,oneOutput)
%       [y,x_,m,k,index,L,U,p] = SIMINIT(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor,oneOutput)
%
% Description:
%       Auxiliary function for the initialization of the simulation functions.       
%
% Input Arguments:
%       -A,B,C,D,E:     state-space matrices
%       -u:             input vector/matrix with dimension Nsample x Ninput
%       -x:             initial state vector for time integration
%       -Ts:            sampling time
%       -Ts_sample:     sampling time for matrix of state-vectors
%       -isDescriptor:  isDescriptor-boolean
%       -oneOutput:     boolean specifying if only y is required as output
%
% Output Arguments:
%       -y:             output vector
%       -x_:            matrix of state vectors
%       -m:             numbers of samples for x
%       -k, index:      initialization of index variables for subsequent
%                       loops
%       -L,U,p:         results of LU factorization
%
% See Also:
%       sim, simForwardEuler, simRK4
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

if oneOutput
    x_      = [];
    m       = inf;
    k       = []; 
    index   = [];
else
    m       = round(Ts_sample/Ts);
    x_      = zeros(length(A),round(size(u,1)/m));
    x_(:,1) = x; %initial state
    k       = 1; %index for state
    index   = 1;
end

% Precompute LU decomposition of E
if isDescriptor
    [L,U,p] = lu(E,'vector');
else
    L = []; U = []; p = [];
end

% Initialize output
y       = zeros(size(C,1),size(u,1));
y(:,1)  = C*x + D*u(1,:)';