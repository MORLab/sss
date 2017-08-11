function [y,x_,tx] = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
% SIMRK4 - Integrates sss model using Runge-Kutta 4
% 
% Syntax:
%       y               = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
%       [y,x_]          = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
%       [y,x_,tx]       = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
%
% Description:
%       Integrates sss model using Runge-Kutta 4
% 
% Input Arguments:       
%       -A,B,C,D,E:     state-space matrices
%       -u:             input vector/matrix with dimension Nsample x Ninput
%       -x:             initial state vector for time integration
%       -Ts:            sampling time
%       -Ts_sample:     sampling time for matrix of state-vectors
%       -isDescriptor:  isDescriptor-boolean
%
% Output Arguments:      
%       -y: output vector
%       -x_: matrix of state vectors
%       -tx: time vector for x_
%
% See Also: 
%       sim, simForwardEuler, simBackwardEuler
%
% References:
%       * *[1] Gear (1971)*, Numerical Initial Value Problems in 
%       Ordinary Differential Equations.
%       * *[2] Shampine (1994)*, Numerical Solution of Ordinary Differential Equations, 
%       Chapman & Hall, New York.
%       * *[3] Shampine and Gordon (1975)*, Computer Solution of Ordinary Differential 
%       Equations: the Initial Value Problem, W. H. Freeman, San Francisco.
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Chair of Thermofluid Dynamics, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Stefan Jaensch, Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  04 Aug 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------ 

%% Initialize variables using common function for all simulation methods
%  Note: using one common simulation function having the methods as nested
%  functions would be much better. Due to historical reasons, the
%  simulation functions not have their present form. Later releases may
%  include some significant restructuring. 

argin = {A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor};

[y,x_,m,k,index,L,U,p] = simInit(argin{:},nargin==1); %common function

%% Run simulation
for i = 2:size(u,1)  
    Bu      = B*u(i-1,:)';    
    % Compute new step
    if ~isDescriptor
        k1  = A*x + Bu;
        k2  = A*(x + Ts/2*k1)   + Bu;
        k3  = A*(x + Ts/2*k2)   + Bu;
        k4  = A*(x + Ts*k3)     + Bu;
    else
        k1  = A*x + Bu;        
        k1  = U\(L\(k1(p,:)));
        k2  = A*(x + Ts/2*k1) + Bu;
        k2  = U\(L\(k2(p,:)));
        k3  = A*(x + Ts/2*k2) + Bu;
        k3  = U\(L\(k3(p,:)));
        k4  = A*(x + Ts*k3) + Bu;
        k4  = U\(L\(k4(p,:)));
    end
    x = (x + Ts/6* (k1 + 2*k2 + 2*k3 + k4));
    
    % Update vectors
    [y,x_,k,index] = simUpdate(y,x_,k,index,x,u,i,m,C,D); %common function
end

tx = (index-1)*Ts;
