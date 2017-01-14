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
% BT tolerance and maximum order for the ROM
tol=1e-6;
max_ord=250;

% ADI tolerance and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.restol = 1e-9;
opts.adi.rctol = 1e-16;
opts.adi.info = 1;
opts.adi.norm = 'fro';

% operations
oper = operatormanager('default');
%%
% Problem data
n0 = 50;  % n0 = number of grid points in either space direction; 
           % n = n0^2 is the problem dimension!
           % (Change n0 to generate problems of different size.)

eqn.A_ = fdm_2d_matrix(n0,'10*x','100*y','0');
eqn.B = fdm_2d_vector(n0,'.1<x<=.3');
eqn.C = fdm_2d_vector(n0,'.7<x<=.9');eqn.C=eqn.C';
n=oper.size(eqn, opts);

eqn.haveE = 0;
%%
%Heuristic Parameters via basic Arnoldi 
opts.adi.shifts.l0=25;
opts.adi.shifts.kp=50;
opts.adi.shifts.km=25;

opts.adi.shifts.b0=ones(n,1);

opts.adi.shifts.p=mess_para(eqn,opts,oper);

disp(opts.adi.shifts.p);
%%
%observability
eqn.type='N';
tic
[ZB,out]=mess_lradi(eqn,opts,oper);
toc
figure(1)
semilogy(out.res);
title('AX + XA^T = -BB^T');
xlabel('number of iterations');
ylabel('normalized residual norm');
pause(1)


disp('size ZB:')
size(ZB)

%%
%controllability
eqn.type = 'T';
tic
[ZC,out]=mess_lradi(eqn,opts,oper);
toc

figure(2)
semilogy(out.res);
title('A^TX + XA = -C^TC');
xlabel('number of iterations');
ylabel('normalized residual norm');
pause(1);

disp('size ZC:')
size(ZC)

%%

[U0,S0,V0] = svd(ZC'*ZB,0);
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
Ar = SC'*(eqn.A_*SB);
Br = SC'*eqn.B;
Cr = eqn.C*SB;

nsample = 200;
w = logspace(-3, 4, nsample);

tr1=zeros(1,nsample); tr2=tr1; err=tr1; relerr=tr1;
fprintf(['Computing TFMs of original and reduced order systems and ' ...
         'MOR errors\n']) 

In=speye(size(eqn.A_,1));
Ir=eye(r);
for k=1:nsample
  fprintf('\r Step %3d / %3d',k,nsample)
  g1 = eqn.C / (1i*w(k)*In - eqn.A_) * eqn.B;
  g2 = Cr / (1i*w(k)*Ir - Ar) * Br;
  err(k) = max(svds(g1-g2));
  tr1(k) = max(svds(g1));
  tr2(k) = max(svds(g2));
  relerr(k)=err(k)/tr1(k);
end
clear I
fprintf('\n\n');

figure(3)
subplot(2,1,1); 
loglog(w, err); 
title('absolute model reduction error')
xlabel('\omega')
ylabel('\sigma_{max}(G(j\omega) - G_r(j\omega))')
subplot(2,1,2); 
loglog(w, relerr);
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
