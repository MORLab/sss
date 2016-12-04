function [Er,Ar,Br,Cr,TL,TR]=mess_balanced_truncation(E,A,B,C,max_order,trunc_tol,info)
% Lyapunov Balanced truncation for descriptor systems with invertible E.
%
%  [Er,Ar,Br,Cr,TL,TR]=mess_balanced_truncation(E,A,B,C,max_order,trunc_tol,info)
%
% INPUTS:
%  E,A,B,C    The mass, system, input and output matrices describing the
%             original system
%  max_ord    maximum reduced order allowed 
%             (optional, defaults to size(A,1))
%  trunc_tol  error tolerance used for the Hankel singular value truncation
%             (optional, defaults to 1e-5)
%  info       verbosity control parameter:
%              0  quiet
%              1  show iteration numbers and residuals
%              >1 plot residual history
%
% OUTPUTS:
% Er,Ar,Br,Cr the reduced order model matrices
% TL, TR      the left and right transformation matrices (optional)
% 

% Author: Jens Saak March 2016

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

% BT tolerance and maximum order for the ROM
if nargin<7
  opts.adi.info=0;
  opts.bt.info=0;
  opts.bt.srm.info=0;
  info=0;
else
  opts.adi.info=info;
  opts.bt.info=info;
  opts.bt.srm.info=info;
end
if (nargin<6 || isempty(trunc_tol)), 
  opts.bt.srm.tol=1e-5; 
else
  opts.bt.srm.tol=trunc_tol; 
end
if ( nargin<5 || isempty(max_order)), 
  opts.bt.srm.max_ord=size(A,1); 
else
  opts.bt.srm.max_ord=max_order; 
end


% ADI tolerance and maximum iteration number
opts.adi.maxiter=100;
opts.adi.restol=min( 1e-9, opts.bt.srm.tol/100 );
opts.adi.rctol=1e-16;
opts.adi.norm='fro';

% operations
oper = operatormanager('default');
%%
% Problem data

if ~issparse(E)||~issparse(A)
    error('MESS:data', 'Both E and A need to be sparse.');
end
if sprank(E) < size(E,1)
    error('MESS:data', 'Only systems with invertible E are supported at the moment')
end

eqn.E_= E;
eqn.A_= A;
eqn.B = B;
eqn.C = C;
eqn.D = [];

n=oper.size(eqn, opts);

if norm(E-speye(n),'inf')==0
    eqn.haveE=0;
else
    eqn.haveE=1;
end


%%
%Heuristic Parameters via basic Arnoldi 
opts.adi.shifts.l0=25;
opts.adi.shifts.kp=50;
opts.adi.shifts.km=25;

opts.adi.shifts.b0=ones(n,1);

opts.adi.shifts.method = 'projection';

opts.adi.shifts.p=mess_para(eqn,opts,oper);

%disp(opts.adi.shifts.p);
%%
%observability
eqn.type='N';
[ZB,out]=mess_lradi(eqn,opts,oper);
if out.niter==opts.adi.maxiter
  warning('MESS:BT',['ADI did not converge for observability Gramian ' ...
                     'factor. Reduction results may be ' ...
                     'inaccurate']);
end

if info>1
  figure(1)
  semilogy(out.res);
  title('AX + XA^T = -BB^T');
  xlabel('number of iterations');
  ylabel('normalized residual norm');
  drawnow
end

if info>0 
  disp('size ZB:')
  size(ZB)
end

%%
%controllability

eqn.type = 'T';
[ZC,out]=mess_lradi(eqn,opts,oper);
if out.niter==opts.adi.maxiter
  warning('MESS:BT',['ADI did not converge for controllability Gramian ' ...
                     'factor. Reduction results may be ' ...
                     'inaccurate']);
end

if info>1
  figure(2)
  semilogy(out.res);
  title('A^TX + XA = -C^TC');
  xlabel('number of iterations');
  ylabel('normalized residual norm');
  drawnow
end
if info>0
  disp('size ZC:')
  size(ZC)
end

%% execute square root method
if nargout == 4
  [Er,Ar,Br,Cr, Dr] = mess_square_root_method(eqn,opts,oper,ZB,ZC);
end
if nargout == 6
  [Er,Ar,Br,Cr,Dr,TL,TR] = mess_square_root_method(eqn,opts,oper,ZB,ZC);
end
if info>2
  mess_sigma_plot(eqn, opts, oper, Er, Ar, Br, Cr, Dr, 1e-6, 1e6, 100)
end