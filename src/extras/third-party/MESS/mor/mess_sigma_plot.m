function mess_sigma_plot(eqn, opts, oper, Er,Ar,Br,Cr,Dr,fmin,fmax,nsample)
% Computation of simple sigma-plots for descriptor systems with invertible
% E.
%
%   mess_sigma_plot(eqn, opts, oper, Er, Ar, Br, Cr, Dr, fmin, fmax, nsample)
%
% INPUTS:
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation 
%                       with A and E
% 
%   Er,Ar,Br,Cr,Dr      reduced order model matrices
%
% fmin, fmax   left and right bounds of the frequency range. They will be
%              interpreted as exponents in the logarithmic range if
%              integers are passed.
% nsample      number of transfer function samples to take in the plot

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

if (floor(fmin)==fmin) && (floor(fmax)==fmax)
    w=logspace(fmin,fmax,nsample);
else
    w=logspace(log10(fmin),log10(fmax),nsample);     
end

tr1=zeros(1,nsample); tr2=tr1; err=tr1; relerr=tr1;

fprintf(['Computing TFMs of original and reduced order systems and ' ...
         'MOR errors\n']) 
%% preprocess shifted solver
[eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);

for k=1:nsample
  if ~mod(k,nsample/10), fprintf('\r Step %3d / %3d',k,nsample); end
  if isfield(eqn,'D')&& ~isempty(eqn.D)
      g1 = eqn.C *oper.sol_ApE(eqn, opts,'N',-1i*w(k),'N',-eqn.B,'N')+eqn.D;
  else
      g1 = eqn.C *oper.sol_ApE(eqn, opts,'N',-1i*w(k),'N',-eqn.B,'N');
  end
  if ~isempty(Dr)
      g2 = Cr / (1i*w(k)*Er - Ar) * Br + Dr;
  else
      g2 = Cr / (1i*w(k)*Er - Ar) * Br;
  end
  err(k) = max(svds(g1-g2));
  tr1(k) = max(svds(g1));
  tr2(k) = max(svds(g2));
  relerr(k)=err(k)/tr1(k);
end
fprintf('\n\n');
%% postprocess shifted solver
[eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper);

figure
subplot(2,1,1); 
loglog(w, err); 
title('absolute model reduction error')
xlabel('\omega')
ylabel('\sigma_{max}(G(j\omega) - G_r(j\omega))')
axis tight
subplot(2,1,2); 
loglog(w, relerr);
title('relative model reduction error')
xlabel('\omega')
ylabel(['\sigma_{max}(G(j\omega) - G_r(j\omega)) / \' ...
        'sigma_{max}(G(j\omega))'])
axis tight
    
figure
loglog(w, tr1)
hold on
loglog(w, tr2, 'r--')
legend({'original system','reduced system'})
xlabel('\omega')
ylabel('\sigma_{max}(G(j\omega))')
title('Transfer functions of original and reduced systems')
axis tight
hold off
