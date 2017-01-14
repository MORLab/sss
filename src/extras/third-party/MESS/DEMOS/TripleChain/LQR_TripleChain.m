%
%%

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
clear all%, close all%, clc

% set operation
oper = operatormanager('so_1');
% Problem data

n1=1000;
alpha=2;
beta=5;
v=5;

[eqn.M_,eqn.D_,eqn.K_]=triplechain_MSD(n1,alpha,beta,v);

s  = size(eqn.K_,1);

eqn.B = [sparse(s,1);ones(size(eqn.K_,1),1)];
eqn.C = eqn.B';

eqn.haveE=1;

%%
% ADI tolerances and maximum iteration number
opts.adi.maxiter = 200;
opts.adi.restol = 1e-10;
%opts.adi.rctol = 1e-16;
opts.adi.rctol = 0;
opts.adi.info = 0;
opts.adi.projection.freq=0;
opts.adi.accumulateK = 1;
opts.adi.accumulateDeltaK = 0;
opts.adi.computeZ = 1;

eqn.type = 'T';
%%
%Heuristic shift parameters via basic Arnoldi 
n=oper.size(eqn, opts);
opts.adi.shifts.l0=25;
opts.adi.shifts.kp=50;
opts.adi.shifts.km=25;

opts.adi.shifts.info=0;
opts.adi.shifts.method = 'heur';
opts.adi.shifts.b0=ones(n,1);
%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 25;
opts.nm.restol = 1e-10;
opts.nm.rctol = 1e-16;
opts.nm.info = 1;
opts.nm.accumulateRes = 1;
opts.nm.linesearch = 1;
opts.nm.inexact = 'quadratic';
opts.nm.tau = 0.1;
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


%%so_2 Transformation Test
%%
clear all, close all

% set operation
oper = operatormanager('so_2');
% Problem data

n1=1000;
alpha=2;
beta=5;
v=5;

[eqn.M_,eqn.D_,eqn.K_]=triplechain_MSD(n1,alpha,beta,v);

s  = size(eqn.K_,1);

eqn.B = [sparse(s,1);eqn.M_\ones(size(eqn.K_,1),1)];
eqn.C = [sparse(s,1);ones(size(eqn.K_,1),1)]';

eqn.haveE=1;

%%

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.restol = 1e-10;
%opts.adi.rctol = 1e-16;
opts.adi.rctol = 0;
opts.adi.info = 0;
opts.adi.projection.freq=0;
opts.adi.accumulateK = 0;
opts.adi.accumulateDeltaK = 1;
opts.adi.computeZ = 0;

eqn.type = 'T';

%%
%Heuristic shift parameters via basic Arnoldi 
n=oper.size(eqn, opts);
opts.adi.shifts.l0=25;
opts.adi.shifts.kp=50;
opts.adi.shifts.km=25;
opts.adi.shifts.method = 'heur';
opts.adi.shifts.info=0;
opts.adi.shifts.b0=ones(n,1);
%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 15;
opts.nm.restol = 1e-10;
opts.nm.rctol = 1e-16;
opts.nm.info = 1;
opts.nm.linesearch = 1;
opts.nm.accumulateRes = 1;
opts.nm.linesearch = 1;
opts.nm.inexact = 'quadratic';
opts.nm.tau = 0.1;
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
