function au_s_i(p)
%
%  Generates the data used in 'au_s'. Data are stored in global
%  variables.
% 
%  NOTE that 'au_m_i' must be called before calling this routine.
% 
%  Calling sequence:
%
%    au_s_i(p)
%
%  Input:
%
%    p         vector containing the ADI shift parameters p(i).
%
%  Remarks:
% 
%    This routine has access to the matrix A, which must be provided
%    by the routine 'au_m_i',
%
%    The real parts of the entries of p must be negative.
%
%
%  LYAPACK 1.6 (Jens Saak, October 2007)

if nargin~=1
  error('Wrong number of input arguments.');
end

if any(real(p)>=0)
  error('Real parts of entries of p must be negative!');
end

l = length(p);

global LP_A LP_LC LP_UC LP_aC LP_oC LP_SC

LP_LC = cell(1,l);
LP_UC = cell(1,l);
LP_aC = cell(1,l);
LP_oC = cell(1,l);
LP_SC = cell(1,l);

if isempty(LP_A)
  error('This routine needs global data which must be generated by calling ''au_m_i'' first.');
end 

I = speye(size(LP_A));

for i = 1:l
%   [LP_LC{i}, LP_UC{i}]= lu(LP_A+p(i)*I);
  [LP_LC{i}, LP_UC{i}, LP_aC{i}, LP_oC{i}, LP_SC{i}] = lu(LP_A+p(i)*I,'vector'); %[L,U,a,o,S]=lu(A-pI),'vector');
  
end
