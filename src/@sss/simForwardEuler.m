function [y,x_,index] = simForwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
% simForwardEuler - Integrates sss model using forward Euler
% 
% Syntax:
%       y = simForwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
%       [y,x_] = simForwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
%       [y,x_,index] = simForwardEuler(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
%
% Description:
%       Integrates sss model using forward Euler
% 
% Input Arguments:       
%       -A,B,C,D,E:     state space matrices
%       -u:             input vector in [Nsample,Ninput]
%       -x:             start vector for time integration
%       -Ts:            Sampling time
%       -Ts_sample:     Sampling time for matrix of state-vectors
%       -isDescriptor:  is descriptor
%
% Output Arguments:      
%       -y: output vector
%       -x_: matrix of state vectors
%       -index: time index for x_
%
% Examples:
%       TODO
%
% See Also: 
%       sss/sim, simBackwardEuler, simRK4
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
% Authors:      Stefan Jaensch, Maria Cruz Varona 
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

y(:,1) = C*x + D*u(1,:)';
ETsA = E+Ts*A; TsB = Ts*B;
if isDescriptor
    [L,U,p] = lu(E,'vector');
end

for i = 2:size(u,1)    
    if ~isDescriptor
        x = ETsA*x + TsB*u(i-1,:)';
    else
        x = ETsA*x + TsB*u(i-1,:)';
        x = U\(L\(x(p,:)));
    end
    y(:,i) = C*x + D*u(i,:)';
    if ~isempty(x_)
        if mod(i,m) == 0
            x_(:,k) = x;
            index = [index i];
            k = k+1;            
        end
    end
end