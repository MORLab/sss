%

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
%
function [y,yr] = impeuler(E,A,B,C,Er,Ar,Br,Cr,tau,tmin,tmax,x0,xr0,alpha)
% Simple implicit Euler implementation for validation of the DAE2
% MESS open loop example via a basic step response computation
%
%  [y,yr] = impeuler(E,A,B,C,Er,Ar,Br,Cr,tau,tmin,tmax,x0,xr0,alpha)
%
% INPUTS:
% E,A,B,C      The original system matrices
% Er,Ar,Br,Cr  The reduced order system matrices
% tau          time step size
% tmin         start time 
% tmax         end time
% x0, x0r      initial states of the full and reduced order models
% alpha        step height for the input step function
%
% OUTPUT:
% y,yr         outputs of the full and reduced systems in [tmin,tmax]
%

% Author: Jens Saak
  [L,U,P,Q] = lu(E-tau*A);
  [Lr,Ur,Pr] = lu(Er-tau*Ar);
 
 ntau=ceil((tmax-tmin)/tau);
 y=zeros(size(C,1),ntau);
 yr=y;
 for i=1:ntau
    if ~mod(i,ceil(ntau/10)), 
      fprintf('\r Implicit Euler step %d / %d',i,ntau); 
    end
    if i<0.1*ntau
        x=Q*(U\(L\(P*(E*x0))));
        xr=Ur\(Lr\(Pr*(Er*xr0)));
    else
        salpha=smoother(alpha,i,ntau);
        x=Q*(U\(L\(P*(E*x0+(salpha*tau)*sum(B,2)))));
        xr=Ur\(Lr\(Pr*(Er*xr0+(salpha*tau)*sum(Br,2))));
    end
    y(:,i)=C*x;
    yr(:,i)=Cr*xr;
    x0=x;
    xr0=xr;
 end
  fprintf('\n\n');
end

function alpha= smoother(alpha,i,ntau)
 if i<0.2*ntau
    alpha=sin((10*pi*(i-0.1*ntau)/ntau)-0.5*pi)*alpha;
 end
end
