function [F,E]=ellip(hk,phi);
%  function [F,E]=ellip(hk,phi);
%  Computes complete and incomplete elliptic integrals F(k,phi) and E(k,phi)
%       Input  : hk  --- Modulus k ( 0 < k < 1 )
%                phi --- Argument 
%       Output : F   --- F(k,phi)
%                E   --- E(k,phi)
%       ==================================================

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

g=0.0;
a0=1.0;
b0=min(1-eps,sqrt(1.0-hk.*hk));
d0=phi;
r=hk.*hk;
if (hk == 1.0&phi == pi/2) ;
  F=1.0e+300;
  E=1.0;
elseif (hk == 1.0);
  F=log((1.0+sin(d0))./cos(d0));
  E=sin(d0);
else;
  fac=1.0;
  for  n=1:40;
    a=(a0+b0)./2.0;
    b=sqrt(a0.*b0);
    c=(a0-b0)./2.0;
    fac=2.0.*fac;
    r=r+fac.*c.*c;
    if (phi ~= pi/2) ;
      d=d0+atan((b0./a0).*tan(d0));
      g=g+c.*sin(d);
      d0=d+pi.*fix(d./pi+.5);
    end;
    a0=a;
    b0=b;
    if (c < 1.0e-15) break; end;
  end;
  ck=pi./(2.0.*a);
  ce=pi.*(2.0-r)./(4.0.*a);
  if (phi == pi/2) ;
    F=ck;
    E=ce;
  else
    F=d0./(fac.*a);
    E=F.*ce./ck+g;
  end;
end;
return;
