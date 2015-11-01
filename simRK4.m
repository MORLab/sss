function [y,x_,index] = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
% Integrates sss model using Runge-Kutta 4
% ------------------------------------------------------------------
% [y,x_,index] = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
% Inputs:       * A,B,C,D,E: state space matrices
%               * u: input vector in [Nsample,Ninput]
%               * x: start vector for time integration
%               * Ts: Sampling time
%               * Ts_sample: Sampling time for matrix of state-vectors
%               * isDescriptor: is descriptor
% Outputs:      * y: output vector
%               * X: matrix of state vectors [Optional]
%               * index: time index for X  [Optional]
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Stefan Jaensch (jaensch@tfd.mw.tum.de)
% Last Change:
% ------------------------------------------------------------------
%
% see also: sss/sim, forwardEuler, backwardEuler


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