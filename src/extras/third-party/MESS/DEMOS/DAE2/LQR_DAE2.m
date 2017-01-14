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
oper = operatormanager('dae_2');

%%
% Problem data
 problem='stokes'
% problem='NSE'
% lvl=1 
% re=300

switch problem
  case 'stokes'
    nin = 5;
    nout = 5;
    nx = 10;
    ny = 10;
    [eqn.E_,eqn.A_,eqn.B,eqn.C]=stokes_ind2(nin,nout,nx,ny);
    eqn.haveE=1;
    st=full(sum(diag(eqn.E_)));
    eqn.st=st;
    eqn.B=eqn.B(1:st,:);
    eqn.C=eqn.C(:,1:st);
  case 'NSE'
    try
      load(sprintf('mat_nse_re_%d',re))
    catch
      error('The files mat_nse_re_300.mat, mat_nse_re_400.mat and mat_nse_re_500.mat ar available for dowload in a separate archive (270MB each). Please fetch them from the MESS download mage and unpack them into the DEMOS/models/NSE folder.')
    end
    eqn.A_=mat.mat_v.fullA{lvl};
    eqn.E_=mat.mat_v.E{lvl};
    eqn.haveE=1;
    eqn.B=mat.mat_v.B{lvl};
    if re>200
      eqn.K0=mat.mat_v.Feed_0{lvl};
      eqn.haveUV=1;
    end
    eqn.C=mat.mat_v.C{lvl};
    eqn.st=mat.mat_mg.nv(lvl);
    st=eqn.st;
end
k=2;

%%

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 300;
opts.adi.restol = 1e-12;
opts.adi.rctol = 1e-16;
opts.adi.info = 1;
opts.adi.projection.freq=0;

eqn.type='T';
%%
n=size(eqn.A_, 1);
opts.adi.shifts.l0=25;
opts.adi.shifts.kp=50;
opts.adi.shifts.km=25;
opts.adi.shifts.method = 'projection';
opts.adi.shifts.l0= 9;
opts.adi.shifts.b0=ones(n,1);
%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 20;
opts.nm.restol = 1e-10;
opts.nm.rctol = 1e-16;
opts.nm.info = 1;
opts.nm.projection.freq=0;
opts.nm.projection.ortho=1;
%opts.nm.projection.meth='care_nwt_fac';
opts.nm.res=struct('maxiter',10,'tol',1e-6,'orth',0);
opts.nm.linesearch = 1;
opts.nm.inexact = 'superlinear';
opts.nm.tau = 0.1;
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
