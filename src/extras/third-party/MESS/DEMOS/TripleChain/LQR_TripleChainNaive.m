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
clear all, close all
%%
% set operation
oper = operatormanager('default');

% Problem data
n1=500;
alpha=2;
beta=5;
v=5;
[M_,D_,K_]=triplechain_MSD(n1,alpha,beta,v);


s = size(K_,1);
eqn.A_ = [sparse(s,s),-K_;-K_,-D_];
eqn.E_ = [-K_,sparse(s,s);sparse(s,s),M_];
eqn.B = [sparse(s,1);ones(size(K_,1),1)];
eqn.C = eqn.B';
nin = 4;
nout = 4;
%eqn.B = [zeros(s, nin);rand(s,nin)];
%eqn.C = [zeros(nout,s),rand(nout,s)];

eqn.haveE=1;

%%
% ADI tolerances and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.restol = 1e-10;
% opts.adi.rctol = 1e-16;
opts.adi.rctol = 0;
opts.adi.info = 1;
opts.adi.projection.freq=0;
opts.adi.accumulateDeltaK = 1;
opts.adi.accumulateK = 0;
opts.adi.computeZ = 1;

eqn.type='T';

%%
%Heuristic shift parameters via basic Arnoldi 
opts.adi.shifts.l0=6;
opts.adi.shifts.kp=50;
opts.adi.shifts.km=25;
n=oper.size(eqn, opts);
opts.adi.shifts.b0=ones(n,1);
opts.adi.shifts.method = 'projection';

%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 15;
opts.nm.restol = 1e-10;
opts.nm.rctol = 1e-16;
opts.nm.info = 1;
opts.nm.accumulateRes = 1;
opts.nm.linesearch = 1;
opts.nm.inexact = 'quadratic';
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
