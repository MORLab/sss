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
clear , close all
%%
% set operation
oper = operatormanager('dae_1');

%% Problem data
load('bips98_606.mat')
% from https://sites.google.com/site/rommes/software
p = find(diag(E));
np = find(diag(E) == 0);
pp = [p;np];
eqn.A_ = A(pp, pp);
eqn.E_ = E(pp, pp);
eqn.B = b(pp, :);
eqn.C = c( : , pp);
eqn.st = length(p);
eqn.haveE = 1;

%%

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 300;
opts.adi.restol = 1e-12;
opts.adi.rctol = 1e-16;
opts.adi.info = 1;
opts.adi.projection.freq=0;

eqn.type='N';
%%
% n=size(eqn.A_, 1);
opts.adi.shifts.l0=25;
opts.adi.shifts.kp=50;
opts.adi.shifts.km=25;
opts.adi.shifts.method = 'projection';
opts.adi.shifts.l0= 9;
% opts.adi.shifts.b0=ones(n,1);
%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 30;
opts.nm.restol = 1e-10;
opts.nm.rctol = 1e-16;
opts.nm.info = 1;
opts.nm.linesearch = 1;
opts.nm.accumulateRes = 1;
opts.nm.norm = 'fro';
%%
tic
[ZB,out]=mess_lrnm(eqn,opts,oper);
toc

figure(1)
out.nm.res
semilogy(out.nm.res);
title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
xlabel('number of iterations');
ylabel('normalized residual norm');
pause(1)

disp('size ZB:')
size(ZB)
