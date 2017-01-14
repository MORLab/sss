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
x0=zeros(size(eqn.A_,1),1);
xr0=zeros(size(Ar,1),1);

alpha=1;
tau=1e-2;
tmin=0;
tmax=50;
T=tmin:tau:tmax;
ntau=floor((tmax-tmin)/tau);
range=1:floor(ntau/500):ntau;
%%
tic
[y,yr] = impeuler_closed(eqn.E_,eqn.A_,eqn.B,eqn.C,eye(size(Ar)),Ar,Br,Cr,K,Kr,tau,tmin,tmax,x0,xr0,alpha);
toc
%%
abserr=abs(y-yr);
relerr=abs(abserr./y);
 
%%
colors=['y','m','c','r','g','b','k'];

figure(20)
hold on
for j=1:size(eqn.C,1)
    plot(T(range),y(j,range),colors(j));
    plot(T(range),yr(j,range),strcat(colors(j),'--'));
end
xlabel('time')
ylabel('magnitude of outputs')
title('step response')
if strcmp(problem,'NSE')
    legend('out1','out1 red','out2','out2 red','out3','out3 red','out4','out4 red','out5','out5 red','out6','out6 red','out7','out7 red','Location','EastOutside')
else
    legend('out1','out1 red','out2','out2 red','out3','out3 red','out4','out4 red','out5','out5 red','Location','EastOutside')
end
hold off
figure(20)
%%
figure(21)

for j=1:size(eqn.C,1)
    semilogy(T(range),abserr(j,range),colors(j));
    if j==1, hold on; end
end
xlabel('time')
ylabel('magnitude')
title('absolute error')
if strcmp(problem,'NSE')
    legend('out1','out2','out3','out4','out5','out6','out7','Location','EastOutside')
else
    legend('out1','out2','out3','out4','out5','Location','EastOutside')
end
hold off

figure(21)

figure(22)

for j=1:size(eqn.C,1)
    semilogy(T(range),relerr(j,range),colors(j));
    if j==1, hold on; end
end
xlabel('time')
ylabel('magnitude')
title('relative error')
if strcmp(problem,'NSE')
    legend('out1','out2','out3','out4','out5','out6','out7','Location','EastOutside')
else
    legend('out1','out2','out3','out4','out5','Location','EastOutside')
end
hold off
figure(22)
