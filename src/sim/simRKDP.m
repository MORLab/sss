function [y,x_,index] = simRKDP(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
% simRKDP - Integrates sss model using Runge-Kutta Dormand�Prince
% 
% Syntax:
%       [y,x_,index] = simRK4(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor)
%
% Description:
%       Integrates sss model using Runge-Kutta Dormand�Prince
%       TODO
%
% Input Arguments:    
%       *Required Input Arguments:*
%           -A,B,C,D,E:     state space matrices
%           -u:             input vector in [Nsample,Ninput]
%           -x:             start vector for time integration
%           -Ts:            Sampling time
%           -Ts_sample:     Sampling time for matrix of state-vectors
%           -isDescriptor:  is descriptor
%
% Output Arguments:      
%       -y: output vector
%       -X: matrix of state vectors [Optional]
% 
% Examples:
%       TODO
%
% See Also: 
%       sss/sim, forwardEuler, backwardEuler* index: time index for X  [Optional]
%
% References:
%       TODO
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
% Authors:      Stefan Jaensch (jaensch@tfd.mw.tum.de)
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  28 Feb 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

warning('This function is work in progress')

Ts0 = Ts;

if nargout == 1
    x_ = [];
else    
    x_ = x;    
    t_save = 0;
    index = [];
end
Bu=B*u(1,:)';
k1 = A*x + Bu;
if isDescriptor
    [L,U,p] = lu(E,'vector');
    k1 = U\(L\(k1(p,:)));
end

y = C*x + D*u(1,:)';
t = 0;
tu = 0:Ts:Ts*(length(u)-1);
i = 2;
while t(end)<=tu(end)
    
    %     Bu = B*u(i-1,:)';
    if ~isDescriptor        
        k2 = A*(x + Ts*k1/5) + Bu;% + B*interp1(tu',u,t(end)+Ts/5)';
        k3 = A*(x + Ts*(3*k1+9*k2)/40) + Bu;% B*interp1(tu,u,t(end)+Ts*3/10)';
        uc = interp1(tu,u,t(end)+Ts)';
        Bu = B*uc;
        k4 = A*(x + Ts*(44/45*k1-56/15*k2+32/9*k3)) + Bu;% B*interp1(tu,u,t(end)+Ts*4/5)';
        k5 = A*(x + Ts*(19372/6561*k1 - 25360/2187*k2 + 64448/6561*k3 -212/729 *k4)) + Bu;% B*interp1(tu,u,t(end)+Ts*8/9)';
        k6 = A*(x + Ts*(9017/3168*k1 - 355/33*k2 + 46732/5247*k3 +49/176*k4 -5103/18656*k5 )) + Bu;% B*interp1(tu,u,t(end)+Ts)';
        x__ = 35/384*k1 + 500/1113*k3 + 125/192*k4 -2187/6784*k5 + 11/84*k6;
        k7 = A*(x + Ts*x__) + Bu;% B*interp1(tu,u,t(end)+Ts)';
        x6thOrder = x+Ts*(5179/57600*k1 + 7571/16695*k3 +  393/640*k4 -92097/339200*k5+ 187/2100*k6+ 1/40 *k7);
        x = x+Ts*x__;
        k1 = k7;
    else
        
        k2 = A*(x + Ts*k1/5) + Bu;% + B*interp1(tu',u,t(end)+Ts/5)';
        k2 = U\(L\(k2(p,:)));
        k3 = A*(x + Ts*(3*k1+9*k2)/40) + Bu;% + B*interp1(tu,u,t(end)+Ts*3/10)';
        k3 = U\(L\(k3(p,:)));
        uc = interp1(tu,u,t(end)+Ts)';
        Bu = B*uc;
        k4 = A*(x + Ts*(44/45*k1-56/15*k2+32/9*k3)) + Bu;% + + B*interp1(tu,u,t(end)+Ts*4/5)';
        k4 = U\(L\(k4(p,:)));
        k5 = A*(x + Ts*(19372/6561*k1 - 25360/2187*k2 + 64448/6561*k3 -212/729 *k4)) + Bu;% B*interp1(tu,u,t(end)+Ts*8/9)';
        k5 = U\(L\(k5(p,:)));
        k6 = A*(x + Ts*(9017/3168*k1 - 355/33*k2 + 46732/5247*k3 +49/176*k4 -5103/18656*k5 )) + Bu;% + B*interp1(tu,u,t(end)+Ts)';
        k6 = U\(L\(k6(p,:)));
        x__ = 35/384*k1 + 500/1113*k3 + 125/192*k4 -2187/6784*k5 + 11/84*k6;
        k7 = A*(x + Ts*x__) + Bu;% +  B*interp1(tu,u,t(end)+Ts)';
        k7 = U\(L\(k7(p,:)));
        x6thOrder = x+Ts*(5179/57600*k1 + 7571/16695*k3 +  393/640*k4 -92097/339200*k5+ 187/2100*k6+ 1/40 *k7);
        x = x+Ts*x__;
        k1 = k7;
    end
    y =[y C*x+D*uc];%interp1(tu,u,t(end)+Ts)'];
    t = [t t(end)+Ts];   
    if ~isempty(x_)
        if t(end) >= (Ts_sample+t_save) || t(end)==tu(end)
            t_save = t(end);
            x_ = [x_ x];
            index = [index length(t)];            
        end
    end
    
    if t(end)==tu(end)
        break;
    end
    %TODO: implement proper criterion for adaptive time stepping
    if norm(x-x6thOrder)/norm(x6thOrder)>1e-6
        Ts = Ts/2;
    else
        Ts = Ts*2;
    end
    i = i+1;
    if t(end)+Ts>tu(end)        
        Ts = tu(end)-t(end);
    end
end
y = full(y);
ti = t(1):Ts0:t(end);
y_ = zeros(size(y,1),length(tu));
for i = 1:size(y,1)
    y_(i,:) = interp1(t,y(i,:),tu);
end
y = y_;
%TODO: define correct output variables
if ~isempty(x_)
    index = interp1(t(index),index,ti,'nearest');
end

