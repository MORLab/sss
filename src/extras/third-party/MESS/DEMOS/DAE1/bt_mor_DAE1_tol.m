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
% Note that we here use the alpha shift suggested by Rommes and coauthors,
% but do not perform the shift back.
p = find(diag(E));
np = find(diag(E) == 0);
pp = [p;np];
eqn.A_ = A(pp, pp)-0.05*E(pp,pp);
eqn.E_ = E(pp, pp);
eqn.B = b(pp, :);
eqn.C = c( : , pp);
eqn.st = length(p);
eqn.haveE = 1;

%%
% BT tolerance and maximum order for the ROM
tol=1e-3;
max_ord=250;

%%
% ADI tolerances and maximum iteration number
opts.adi.maxiter = 300;
opts.adi.restol = 1e-12;
opts.adi.rctol = 1e-16;
opts.adi.info = 1;
opts.adi.projection.freq=0;
opts.adi.norm = 'fro';

%%
opts.adi.shifts.kp = 50;
opts.adi.shifts.km = 35;
opts.adi.shifts.method = 'projection';
opts.adi.shifts.l0= 20;

opts.adi.shifts.p=mess_para(eqn,opts,oper);

disp(opts.adi.shifts.p);

%%
eqn.type='N';
tic
[ZB,out]=mess_lradi(eqn,opts,oper);
toc

figure
disp('normalize residuals')
semilogy(out.res);
title('0= BB^T + AXM^T + MXA^T');
xlabel('number of iterations');
ylabel('normalized residual norm');

disp('size ZB:')
size(ZB)

%%
eqn.type='T';
tic
[ZC,out]=mess_lradi(eqn,opts,oper);
toc

figure
disp('normalize residuals')
semilogy(out.res);
title('0= C^TC + A^TXE + E^TXA');
xlabel('number of iterations');
ylabel('normalized residual norm');
pause(1)

disp('size ZC:')
size(ZC)
%%

[U0,S0,V0] = svd(ZC'*(eqn.E_(1:eqn.st,1:eqn.st)*ZB),0);
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

b1 = SC'*(eqn.A_(1:eqn.st,1:eqn.st))*SB;
b2 = SC'*(eqn.A_(1:eqn.st,eqn.st+1:end));
a1 = eqn.A_(eqn.st+1:end,1:eqn.st)*SB;

Ar =  b1 - b2*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\a1);%+a1*eye(r,r);
Br = SC'*eqn.B(1:eqn.st,:) - b2*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\eqn.B(eqn.st+1:end,:));
Cr = eqn.C(:,1:eqn.st)*SB - eqn.C(:,eqn.st+1:end)*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\a1);
Dr = -eqn.C(:,eqn.st+1:end)*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\eqn.B(eqn.st+1:end,:));
%%
nsample = 200;
w = logspace(-3, 4, nsample);

tr1=zeros(1,nsample); tr2=tr1; err=tr1; relerr=tr1;
fprintf(['Computing TFMs of original and reduced order systems and ' ...
         'MOR errors\n']) 

Ir=eye(r);
for k=1:nsample
  if mod(k,10)==0, fprintf('\r Step %3d / %3d',k,nsample), end
  g1 = eqn.C / (1i*w(k)*eqn.E_ - eqn.A_) * eqn.B;
  g2 = Cr / (1i*w(k)*Ir - Ar) * Br + Dr;
  err(k) = max(svds(g1-g2));
  tr1(k) = max(svds(g1));
  tr2(k) = max(svds(g2));
  relerr(k)=err(k)/tr1(k);
end
clear I
fprintf('\n\n');

figure
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


figure
loglog(w, tr1)
hold on
loglog(w, tr2, 'r--')
legend({'original system','reduced system'})
xlabel('\omega')
ylabel('\sigma_{max}(G(j\omega))')
title('Transfer functions of original and reduced systems')
hold off


figure
semilogy(diag(S0));
title('Computed Hankel singular values');
xlabel('index');
ylabel('magnitude');

