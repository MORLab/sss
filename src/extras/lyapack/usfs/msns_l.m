function Y = msns_l(tr,X)
%
%  Solves linear systems with the symmetric, negative definite matrix A, 
%  i.e., Y = inv(A)*X.
%
%  A is given implicitely as A = inv(MU')*N*inv(MU). MU and the Cholesky 
%  factor of N are provided as global data. These data must be generated 
%  by calling 'msns_l_i' before calling this routine!
%  
%  Calling sequence:
%
%    Y = msns_l(tr,X)
%
%  Input:
%
%    tr        is not referenced;
%    X         matrix of proper size.
%
%  Output:
%
%    Y         the solution matrix. 
%
% 
%   LYAPACK 1.0 (Thilo Penzl, May 1999)

if nargin~=2
  error('Wrong number of input arguments.');
end

global LP_NU LP_MU

if isempty(LP_NU) || isempty(LP_MU)
  error('This routine needs global data which must be generated by calling ''msns_l_i'' first.');
end 

Y = -LP_MU*(LP_NU\(LP_NU'\(LP_MU'*X)));      % Note the minus!

