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
oper = operatormanager('dae3_so');

%% load problem data
load('g600.mat');
% load('g6000.mat');
eqn.M_=M;
eqn.D_=D;
eqn.K_=K;
eqn.G_=G;
eqn.haveE=1;
eqn.alpha = -0.02;
nv = size(eqn.M_,1);
np = size(eqn.G_,1);
eqn.B = B(1:2*nv,:);
eqn.C = C(:,1:2*nv);
eqn.type = 'T';

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

eqn.type = 'N';
tic
[Z,out]=mess_lradi(eqn,opts,oper);
toc

opts.adi.info = 0;

%% lrnm

eqn.type = 'T';
tic
[ZB,out]=mess_lrnm(eqn,opts,oper);
toc
