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

k=1;
[eqn.E_,eqn.A_,eqn.B,eqn.C]=getrail(k);
eqn.haveE=1;
%%

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.restol = 1e-12;
opts.adi.rctol = 1e-16;
opts.adi.info = 1;

eqn.type = 'T';


%%
%Heuristic shift parameters via basic Arnoldi 
n=oper.size(eqn, opts);
opts.adi.shifts.kp=50;
opts.adi.shifts.km=25;
opts.adi.shifts.l0 = 6;

opts.adi.shifts.b0=ones(n,1);
opts.adi.shifts.method = 'projection';
opts.adi.shifts.p=mess_para(eqn,opts,oper);
opts.adi.norm = 'fro';

%%
fprintf('########################\n')
fprintf('# ADI without LDL^T:\n')
fprintf('########################\n')
tic
[Z,out]=mess_lradi(eqn,opts,oper);
toc

figure(1)
semilogy(out.res);
xlabel('number of iterations');
ylabel('normalized residual norm');
pause(1)

disp('size ZB:')
size(Z)

%% Set LDL fields
opts.adi.LDL_T = 1;
eqn.S = diag([4,4,9,9,16,16]);
eqn.G = eqn.C' * diag([4,4,9,9,16,16].^(-0.5));

%%
fprintf('########################\n')
fprintf('# ADI with LDL^T:\n')
fprintf('########################\n')
opts.adi.shifts.p=mess_para(eqn,opts,oper);
tic
[Z1,out1]=mess_lradi(eqn,opts,oper);
toc

figure(2)
semilogy(out1.res);
xlabel('number of iterations');
ylabel('normalized residual norm');
pause(1)

disp('size ZB1:')
size(Z1)

%% Error of Lypunov solution
err = norm(Z * Z' - Z1 * kron(diag(out1.D), eqn.S) * Z1') / norm(Z * Z');
fprintf('Relative difference between solution with and without LDL^T: \t %g\n', err)
