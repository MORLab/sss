%
%% clear all close all, clc

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
clear all, close all, clc

%% set operation
oper = operatormanager('dae2_so');

%% load problem data
% generate problem
nv = 500; np=10; nin=2; nout=3; p=0.2;
alpha = 0.1; beta = 0.1; v = 5;

[M, D, K]=triplechain_MSD(nv,alpha,beta,v);
nv=size(M,1);
B = rand(nv,nin);
C = rand(nout,nv);
G = zeros(nv,np);
while rank(full(G))~=np, G = sprand(np,nv,p);end 
eqn.M_=M;
eqn.D_=-D;
eqn.K_=-K;
eqn.G_=G;
eqn.haveE=1;
eqn.alpha = -0.02;
eqn.B = [zeros(nv, nin);rand(nv,nin)];
eqn.C = [zeros(nout,nv),rand(nout,nv)];
eqn.type = 'N';

%% condition numbers of generated input data

A = [zeros(nv,nv),eye(nv,nv),zeros(nv,np); ...
     K, D, G';... 
     zeros(np,nv), G, zeros(np,np)];

E = [eye(nv,nv), zeros(nv,nv), zeros(nv,np);...
     zeros(nv,nv), M, zeros(nv,np); ...
     zeros(np,2*nv+np)];

As = full([zeros(nv,nv),eye(nv,nv),zeros(nv,np); ...
     K, D, G';... 
     G,zeros(np,nv), zeros(np,np)]);

Es = full([eye(nv,nv), zeros(nv,nv), zeros(nv,np);...
     zeros(nv,nv), M, eqn.alpha*G';...
     zeros(np,nv),eqn.alpha*G,zeros(np,np)]);

fprintf('cond(M)=%e\n',condest(eqn.M_));
fprintf('cond(D)=%e\n',condest(eqn.D_));
fprintf('cond(K)=%e\n',condest(eqn.K_));
fprintf('cond(A)=%e\n',condest(As));

%% options
% Adi options
opts.adi.maxiter = 1000;
opts.adi.restol = 1e-15;
opts.adi.rctol = 1e-16;
opts.adi.info = 0;
opts.adi.projection.freq=0; 
opts.adi.shifts.l0=25;
opts.adi.shifts.kp=50;
opts.adi.shifts.km=25;
opts.adi.shifts.b0=ones(2*nv+np,1);
opts.adi.norm = 'fro';
opts.adi.shifts.method = 'projection';
opts.adi.shifts.l0=6;

% newton options and maximum iteration number
opts.nm.maxiter = 20;
opts.nm.restol = 1e-10;
opts.nm.rctol = 1e-16;
opts.nm.info = 1;
opts.nm.accumulateRes = 0;
opts.nm.linesearch = 1;
opts.nm.norm = 'fro';
opts.nm.projection.freq=0;
opts.nm.projection.ortho=1;
opts.nm.res=struct('maxiter',10,'tol',1e-6,'orth',0);

%% lradi

eqn.type = 'T';
opts.adi.info = 1;
tic
opts.adi.shifts.p = mess_para(eqn, opts, oper);

[Z,out]=mess_lradi(eqn,opts,oper);
toc
opts.adi.info = 0;


%% lrnm

eqn.type = 'T';
tic
[ZB,out]=mess_lrnm(eqn,opts,oper);
toc

%non transpose equation
%tranpose your data  
eqn.M_ = eqn.M_';
eqn.D_ = eqn.D_';
eqn.K_ = eqn.K_';
eqn.B = eqn.C';
eqn.C = eqn.B';
tic
[ZB,out]=mess_lrnm(eqn,opts,oper);
toc 
