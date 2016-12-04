function p=mess_wachspress(a,b,alpha,TOL)
%
% function p=mess_wachspress(a,b,alpha,TOL) 
%  
% calculates the optimal ADI shiftparameters (for equations where
% Matrix A is stable and symmetric) according to Jing-Rebecca Li
% and Jakob Whites "Low Rank Solution of Lyapunov equation" (which
% gives an overview of Wachspress`s method form e.g. "The ADI model Problem" 
%
% p   is the array of shift parameters
%  
% a   is assumed to be the absolute value of the smallest magnitude
%     eigenvalue of the Systemmatrix A 
%
% b   is assumed to be the absolute value of the largest magnitude eigenvalue
%
% alpha is the arctan of the maximum of abs(imag(lamba))/abs(real(lambda))
%       for all lambda in the spectrum of A
%
% TOL is the epsilon1 of the above paper. The smaller the better
%     the approximation but the more parameters are calculated 
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
if ~isnumeric(a) || (length(a) ~= 1)
    error('MESS:error_arguments','a has to be a numeric value')
end
if ~isnumeric(b) || (length(b) ~= 1)
    error('MESS:error_arguments','b has to be a numeric value')
end
if ~isnumeric(alpha) || (length(alpha) ~= 1)
    error('MESS:error_arguments','alpha has to be a numeric value')
end
if ~isnumeric(TOL) || (length(TOL) ~= 1)
    error('MESS:error_arguments','TOL has to be a numeric value')
end

  if (alpha==0)
    kprime=a/b;
  else
    c2 = 2/(1+(a/b+b/a)/2);
    m = 2*cos(alpha)*cos(alpha)/c2 -1;
    if (m<1)
      error(['Shift parameters would be complex, function not aplicable, ' ...
             'aborting!']); 
      
      %
      % FIX ME: if m<1 parameter become complex! switch back to the
      % heuristics by Thilo or complex parameters suggested by
      % Wachspress.
      % 
      % ALA2006 -> also test switching to Krylov projection based
      % method (see V.Simoncini)
      %
    end
    kprime = 1/(m+sqrt(m^2-1));
  end
  k=min(1-eps,sqrt(1-kprime^2));
  %this is a workaround for the case
  %k=1 that works for our model reduction problems 
  % (not really nice but it works great for now).
			 
			
%TODO: check the computation of k, kprime to avoid roundoff errors
%and probably replace the hack above. 
				
  [K,E]=ellip(k,pi/2);
  if (alpha==0)
    [v,E]=ellip(kprime,pi/2);
  else
    [v,E]=ellip(kprime,asin(sqrt(a/(b*kprime))));
  end
  J=ceil(K/(2*v*pi)*log(4/TOL));

  p=ones(J,1);
  for i=1:J
    p(i)=-sqrt(a*b/kprime)*dn((i-0.5)*K/J,k); 
                               %here we have the choice to take the
                               %matlab function ellipj or our own
                               %one dn. the later can be proted to
                               %FORTRAN or C Code very easily
  end

  
