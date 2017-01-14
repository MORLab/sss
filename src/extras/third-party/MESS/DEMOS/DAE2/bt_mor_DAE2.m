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

% Computes a standard ROM by solving the generalized lyapunov equations and
% regaining the observability Gramian for the equivalent standard state
% space system from that for the generalized system by multiplication with
% M (see line 112).

% BT tolerance and maximum order for the ROM
tol=1e-5;
max_ord=250;

% ADI tolerance and maximum iteration number
opts.adi.maxiter = 350;
opts.adi.restol = sqrt(eps);
opts.adi.rctol = 1e-16;
opts.adi.info = 1;
opts.adi.shifts.info = 1;
opts.adi.norm = 'fro';

oper = operatormanager('dae_2');
%%
% Problem data
problem='stokes'

%problem='NSE'
%lvl=1
%re=300

switch problem
  case 'stokes'
    nin = 5;
    nout = 5;
    nx = 10;
    ny = 10;
    [eqn.E_,eqn.A_,eqn.B,eqn.C]=stokes_ind2(nin,nout,nx,ny);
    n=size(eqn.E_,1);
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
    eqn.B=mat.mat_v.B{lvl};
    eqn.C=mat.mat_v.C{lvl};
    eqn.st=mat.mat_mg.nv(lvl);
    st=eqn.st;
    eqn.haveE=1;
    n=size(eqn.E_,1);
end
%%
eqn.type='N';
eqn.G=eqn.B;
if strcmp(problem,'NSE') && (re>200)
    eqn.V=mat.mat_v.Feed_1{lvl};
    eqn.U = eqn.C';
    eqn.haveUV=1;
end


opts.adi.shifts.l0=6;
opts.adi.shifts.kp=40;
opts.adi.shifts.km=40;
opts.adi.shifts.method='projection';

opts.adi.shifts.b0=ones(size(eqn.A_,1),1);

opts.adi.shifts.p=mess_para(eqn,opts,oper);

tic
[ZB,out]=mess_lradi(eqn,opts,oper);
toc

figure(1)
semilogy(out.res);
title('AXM^T + MXA^T = -BB^T');
xlabel('number of iterations');
ylabel('normalized residual norm');
pause(1)

disp('size ZB:')
size(ZB)

%%
eqn.type = 'T';
eqn.G=eqn.C';
if strcmp(problem,'NSE') && (re>200)
    if (re == 500) && (lvl == 1)
        eqn.V=mat.mat_v.Feed_0{lvl}';
    else
        eqn.V=mat.mat_v.Feed_0{lvl};
    end
    eqn.U = eqn.B;
    eqn.haveUV=1;
end


opts.adi.shifts.l0=6;
opts.adi.shifts.kp=40;
opts.adi.shifts.km=40;
opts.adi.shifts.method='projection';

opts.adi.shifts.b0=ones(size(eqn.A_,1),1);

opts.adi.shifts.p=mess_para(eqn,opts,oper);

tic
[ZC,out]=mess_lradi(eqn,opts,oper);
toc

figure(2)
semilogy(out.res);
title('A^TXM + M^TXA = -C^TC');
xlabel('number of iterations');
ylabel('normalized residual norm');
pause(1);

disp('size ZC:')
size(ZC)

%%
tic 
[U0,S0,V0] = svd(ZC'*(eqn.E_(1:st,1:st)*ZB),0);
s0=diag(S0);
ks=length(s0);
k=ks;
while (sum(s0(k-1:ks))<tol/2)&&(k>2)
  k=k-1; 
end
k0=k

r= min([max_ord k0]);
fprintf(1,'reduced system order: %d\n\n',r);

sigma_0=diag(S0);
sigma_r=diag(S0(1:r,1:r));

VB = ZB*V0(:,1:r);
VC = ZC*U0(:,1:r);

SB = VB*diag(ones(r,1)./sqrt(sigma_r));
SC = VC*diag(ones(r,1)./sqrt(sigma_r));
Ar = SC'*(eqn.A_(1:st,1:st)*SB);
Br = SC'*eqn.B(1:st,:);
Cr = eqn.C(:,1:st)*SB;
toc
%%
tic
nsample = 200;
w = logspace(-3, 4, nsample);

tr1=zeros(1,nsample); tr2=tr1; err=tr1; relerr=tr1;
fprintf(['Computing TFMs of original and reduced order systems and ' ...
         'MOR errors\n']) 

eqn.B=[eqn.B; sparse(n-st,size(eqn.B,2))];
eqn.C=[eqn.C, sparse(size(eqn.C,1),n-st)];
for k=1:nsample
  if ~mod(k,nsample/10), fprintf('\r Step %3d / %3d',k, nsample); end
  g1 = eqn.C / (1i*w(k)*eqn.E_ - eqn.A_) * eqn.B;
  g2 = Cr / (1i*w(k)*eye(r) - Ar) * Br;
  err(k) = max(svds(g1-g2));
  tr1(k) = max(svds(g1));
  tr2(k) = max(svds(g2));
  relerr(k)=err(k)/tr1(k);
end
fprintf('\n\n');
toc
%%
figure(3)
subplot(2,1,1); 
loglog(w, err); 
hold on 
loglog(w,tol*ones(size(w)),'r--')
hold off
title('absolute model reduction error')
xlabel('\omega')
ylabel('\sigma_{max}(G(j\omega) - G_r(j\omega))')
subplot(2,1,2); 
loglog(w, relerr);
hold on 
loglog(w,tol*ones(size(w)),'r--')
hold off
title('relative model reduction error')
xlabel('\omega')
ylabel(['\sigma_{max}(G(j\omega) - G_r(j\omega)) / \' ...
        'sigma_{max}(G(j\omega))'])

figure(4)
loglog(w, tr1)
hold on
loglog(w, tr2, 'r--')
legend({'original system','reduced system'})
xlabel('\omega')
ylabel('\sigma_{max}(G(j\omega))')
title('Transfer functions of original and reduced systems')
hold off

figure(5)
semilogy(diag(S0));
title('Computed Hankel singular values');
xlabel('index');
ylabel('magnitude');

%%
fprintf(['\nComputing open loop step response of original and reduced order ' ...
         'systems and time domain MOR errors\n']) 
open_step
%%
fprintf('\nComputing ROM based feedback\n') 

[~,~,Kr]=care(Ar,Br,Cr'*Cr);
K=[Kr*SC'*eqn.E_(1:st,1:st),zeros(size(Kr,1),n-st)];
%%
fprintf(['\nComputing closed loop step response of original and reduced order ' ...
         'systems and time domain MOR errors\n']) 
closed_step
