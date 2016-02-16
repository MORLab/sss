function Y = as_s(tr,X,i)
%
%  Solves shifted linear systems with the real, symmetric, negative definite 
%  matrix A, i.e., Y = inv(A+p(i)*I)*X.
%
%  The Cholesky factor of -A-p(i)*I is provided as global data. This data 
%  must be generated by calling 'as_s_i' before calling this routine!
%
%  Calling sequence:
%
%    Y = as_s(tr,X,i)
%
%  Input:
%
%    tr        is not referenced;
%    X         a matrix of proper size;
%    i         the index of the shift parameter.
%
%  Output:
%
%    Y         the resulting solution matrix.
%  
%
%   LYAPACK 1.6 (Jens Saak, Octber 2007)

if nargin~=3
  error('Wrong number of input arguments.');
end

global LP_UC

is_init = length(LP_UC{i});
if ~is_init
  error('This routine needs global data which must be generated by calling ''as_s_i'' first.');
end 

Y = -LP_UC{i}\(LP_UC{i}'\X);      % Note the minus!




