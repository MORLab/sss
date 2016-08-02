function [eqn, erg] = init_lu(eqn, opts,flag1,flag2)
%
%  Preprocessing of the system
%    .
%    x  =  A x + B u
%    y  =  C x,
%
%  where A is SPARSE.
%
%  The preprocessing consists of a permutation of the state 
%
%    x <-- P * x
%
%  with a permutation matrix P for bandwidth reduction, which 
%  results in "overwriting" the system matrices as
%
%    A <-- P * A * P',  B <-- P * B,  C <-- C * P'.
%
%  The bandwidth of the reordered matrix A is often much smaller than that
%  of the original matrix. 
%
%  Note that this preprocessing does not affect the input-output
%  mapping of the dynamical system.
%
%  This routine can also be applied when there is no underlying dynamical
%  system. For example, this is the case when only the Lyapunov equation
%
%    A*X + X*A' = - B*B'   ( or A'*X + X*A = - C'*C )
%   
%  needs to be solved. Here, B (or C) can be omitted; see (1) (or (2)).
%
%  Calling sequence:
%
%    [A,B,C,prm,iprm] = au_pre(A,B,C)
%    [A,dummy,C,prm,iprm] = au_pre(A,[],C)                       (1)
%    [A,B,dummy,prm,iprm] = au_pre(A,B,[])                       (2)
%
%  Input:
%
%    A         n-x-n system matrix; 
%    B         n-x-m system matrix;
%    C         q-x-n system matrix.
%
%  Output:
%
%    A, B, C   permuted system matrices;
%    prm       the permutation that has been used;
%    iprm      the inverse permutation (needed to re-reorder certain data
%              in postprocessing);
%    dummy     a dummy output argument (dummy = [] is returned).
%
%
%  LYAPACK 1.0 (Thilo Penzl, May 1999)

% Input data not completely checked!

%start checking
na = nargin;
if(na<=2)
    error('MESS:control_data','Number of input Arguments are at least 3');

%erg = init_default(eqn, flag1);    
elseif(na==3)
    switch flag1
        case {'A','a'}
            [eqn,erg] = checkA(eqn);
        case {'E','e'}
            [eqn,erg] = checkE(eqn);
        otherwise
            error('MESS:control_data','flag1 has to be ''A_'' or ''E_''');
    end
    
%erg = init_default(eqn,flag1,flag2);
elseif(na==4)
    switch flag1
        case {'A','a'}
            [eqn, erg] = checkA(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA] = checkA(eqn);
                    erg = erg && ergA;
                case {'E','e'}
                    [eqn, ergE]= checkE(eqn);
                    erg = erg && ergE; 
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        case {'E','e'}
            [eqn, erg] = checkE(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA] = checkA(eqn);
                    erg = erg && ergA;
                case {'E','e'}
                    [eqn, ergE]= checkE(eqn);
                    erg = erg && ergE; 
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end 
end
end

%checkdata for A_
function [eqn, erg] = checkA(eqn)
erg = isfield(eqn,'A_');
if(erg)
    erg = isnumeric(eqn.A_);
end
erg=erg&&(size(eqn.A_,1)==size(eqn.A_,2));
end

%checkdata for E_
function [eqn, erg] = checkE(eqn)
if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
if ~eqn.haveE
  erg = 1;
  eqn.E_= speye(size(eqn.A_,1)); %make sure we have an identity for
                                 %computations in ApE functions
else
  erg = isfield(eqn,'E_');
  if(erg)
    erg = isnumeric(eqn.E_);
  end
  erg=erg&&(size(eqn.E_,1)==size(eqn.E_,2));
end
end


