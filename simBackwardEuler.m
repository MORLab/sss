function [y,x_,index] = simBackwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample)
% Integrates sss model using backward Euler
% ------------------------------------------------------------------
% [y,x_,index] = simBackwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
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
% see also: sss/sim, forwardEuler, RK4

y = zeros(size(C,1),size(u,1));
if nargout == 1
    x_ = [];
else
    m = round(Ts_sample/Ts);
    x_ = zeros(length(A),round(size(u,1)/m));    
    k = 1;
    index = [];
end

y(:,1) = C*x + D*u(1,:)';
ETsA = E-Ts*A; TsB = Ts*B;
[L,U,p] = lu(ETsA,'vector');

for i = 2:size(u,1)
%     x = ETsA\(E*x + TsB*u(i,:)');
    g = E*x + TsB*u(i,:)';
    x = U\(L\(g(p,:)));

    y(:,i) = C*x + D*u(i,:)';
    if ~isempty(x_)
        if mod(i,m) == 0
            x_(:,k) = x;
            index = [index i];
            k = k+1;            
        end
    end
end