function Y = munu_m(tr,X)
%
%  Evaluates matrix-matrix products with the real matrix A or its 
%  transposed A':
%
%  for tr = 'N':
%
%    Y = A*X,
%
%  for tr = 'T':
%
%    Y = A'*X.
%
%  A is given implicitely as A = inv(ML')*N*inv(MU). ML, MU and N are
%  provided as global data. These data must be generated by calling
%  'munu_m_i' before calling this routine!
%  
%  If called without input parameters, this routine returns the order of
%  the matrix A.
%
%  Calling sequence:
%
%    Y = munu_m(tr,X)
%    n = munu_m
%
%  Input:
%
%    tr        (= 'N' or 'T') determines whether products with A or A' 
%              should be computed;
%    X         a matrix of proper size.
%
%  Output:
%
%    Y         the resulting product;
%    n         the order of the matrix A.
%
%
%  LYAPACK 1.0 (Thilo Penzl, September 1999)

ni = nargin;

if ni~=2 && ni~=0
  error('Wrong number of input arguments.');
end

global LP_ML LP_MU LP_N

if isempty(LP_ML) || isempty(LP_MU) || isempty(LP_N)
  error('This routine needs global data which must be generated by calling ''munu_m_i'' first.');
end 

if ni==0
  Y = size(LP_N,1);  
else
  if tr=='N'
    Y = LP_ML\(LP_N*(LP_MU\X));
  elseif tr=='T'
    Y = LP_MU.'\(LP_N.'*(LP_ML.'\X));
  else
    error('tp must be either ''N'' or ''T''.');
  end
end

