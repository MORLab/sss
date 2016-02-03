function [y,x_,tx] = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
% SIMRK4 - Integrates sss model using Runge-Kutta 4
% 
% Syntax:
%       y = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
%       [y,x_] = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
%       [y,x_,index] = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
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
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Stefan Jaensch 
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------ 

y = zeros(size(C,1),size(u,1));
if nargout == 1
    x_ = [];
else
    m = round(Ts_sample/Ts);
    x_ = zeros(length(A),round(size(u,1)/m));    
    k = 1;
    index = [];
end
if isDescriptor
    [L,U,p] = lu(E,'vector');
end

y(:,1) = C*x + D*u(1,:)';
for i = 2:size(u,1)  
    
    Bu = B*u(i-1,:)';    
    if ~isDescriptor
        k1 = A*x + Bu;
        k2 = A*(x + Ts/2*k1) + Bu;
        k3 = A*(x + Ts/2*k2) + Bu;
        k4 = A*(x + Ts*k3) + Bu;
    else
        k1 = A*x + Bu;        
        k1 = U\(L\(k1(p,:)));
        k2 = A*(x + Ts/2*k1) + Bu;
        k2 = U\(L\(k2(p,:)));
        k3 = A*(x + Ts/2*k2) + Bu;
        k3 = U\(L\(k3(p,:)));
        k4 = A*(x + Ts*k3) + Bu;
        k4 = U\(L\(k4(p,:)));
    end
    x = (x + Ts/6* (k1 + 2*k2 + 2*k3 + k4));
    y(:,i) = C*x + D*u(i,:)';
    if ~isempty(x_)
        if mod(i,m) == 0
            x_(:,k) = x;
            index = [index i];
            k = k+1;            
        end
    end
end
if m==inf
    x_ = x;
    index = size(u,1);
end
if nargout>1
    tx = index*Ts;
end