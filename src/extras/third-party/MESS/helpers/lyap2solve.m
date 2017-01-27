function X = lyap2solve(A,B)
% Solve Lyapunov equation AX+XA^T+B^T=0
% 
%  Solve Lyapunov equation AX+XA^T+B^T=0
%  via Zhou and Sorensen 2-solve  method  
%    
% Input:
%  A         Matrix from Lyapunov equation
%  B         Matrix from Lyapunov equation
%
% Output:
%  X         Solution of Lyapunov equation

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

m = size(A,1); n = size(B,2);
[Q,R]=schur(A);
idx = m:-1:1; Q2=Q(:,idx); R2=R(idx,idx)';
B = Q'*B*Q2;

Rsq = R*R; I=speye(m); 
j=1;

X=zeros(size(A));

while (j < n+1)
 
  if j==n || (j<n  && abs(R2(j+1,j))<10*eps*max( abs(R2(j,j)), abs(R2(j+1,j+1))) )
 
      if (j>1), b = -B(:,j) - X(:,1:j-1)*R2(1:j-1,j);else
              b = -B(:,j) ;end
      X(:,j) = (R+R2(j,j)*I)\b;
      j = j +1;
  else

      r11 = R2(j,j);  r12 = R2(j,j+1);
      r21 = R2(j+1,j); r22 = R2(j+1,j+1);

      if (j>1), b = -B(:,j:j+1) - X(:,1:j-1)*R2(1:j-1,j:j+1); else
                b = -B(:,j:j+1); end
      b = [R*b(:,1)+r22*b(:,1)-r21*b(:,2), R*b(:,2)+r11*b(:,2)-r12*b(:,1)];
      X(:,j:j+1) = ( Rsq+(r11+r22)*R + (r11*r22-r12*r21)*I)\b;
      j = j + 2;
  end
end

X = Q*X*Q2';


