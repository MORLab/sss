function [E,A,B,C,M,D,K,G,Pl,Pr] = msd_ind3(g,mas,k1,k2,d1,d2)
%
% Damped mass-spring system with a holonomic constraint
%
%       E*x'(t) = A*x(t) + B*u(t),
%          y(t) = C*x(t),
%
% with 
%  E = [ I 0 0 ], A = [ 0  I  0  ], B = [0 ], C =[ C1, C2, C3 ],
%      [ 0 M 0 ]      [ K  D -G' ]      [B2]
%      [ 0 0 0 ]      [ G  0  0  ]      [0 ]
% where
%    M is the mass matrix (M is diagonal, positive definite)
%    K is the stiffness matrix (K is three diagonal)
%    D is the damping matrix (D is three diagonal)
%    G is the constraint matrix of full row rank (=> s*E-A is of index 3)
%     ( here G = [ 1, 0, ..., 0, -1 ] )
%
% INPUT:
%      g  = the number of masses ( g>1 )
%    mas  = the g-by-1 array of masses
%     k1  = the (g-1)-by-1 array of spring constants of the first type
%     k2  = the g-by-1 array of spring constants of the second type
%     d1  = the (g-1)-by-1 array of spring constants of the first type
%     d2  = the g-by-1 array of spring constants of the second type 
%
% Note: 1) if mas=[], k1=[], k2=[], d1=[], d2=[], then the matrices 
%          M,K and D are as in [1]. 
%       2) Elements of mas, k1, k2, d1, d2 should be positive
%
% OUTPUT:
%      E   real n-by-n sparse matrix with 
%           n = (nx+1)*ny+nx*(ny+1)+(nx+1)*ny+nx*(ny+1)
%      A   real n-by-n sparse matrix 
%      B   real n-by-m sparse matrix, m is the number of inpits 
%           ( here B = e_{g+1} )
%      C   real p-by-n sparse matrix, p is the number of outputs
%           ( here C = [ e_1, e_2, e_{q-1} ]' )
%     Pr   the spectral projector onto the right deflating subspace 
%          corresponding to the finite eigenvaluesof s*E-A (sparse)
%     Pl   the spectral projector onto the left deflating subspace 
%          corresponding to the finite eigenvaluesof s*E-A (sparse)
%
% DESCIPTION:
%     The i-th mass of weight m_i is connected to the (i+1)-st mass
%     by a spring and a damper with constants k1_i and d1_i, respectively,
%     and also to the ground by a spring and a damper with constants 
%     k2_i and d2_i, respectively. Additionally, the first mass is 
%     connected to the last one by a rigid bar and it is influenced by 
%     the control u(t). The position of the 1-st, 2-nd and (g-1)-st masses
%     is measured.
%
%     |-------------------------------------------------------------------------------|
%     |                                                                               |
%     |      k1_1                  k1_i               k1_i+1                k1_g-1    |
%     |  |-\/\/\/\/-|         |-\/\/\/\/-|         |-\/\/\/\/-|         |-\/\/\/\/-|  |
% m_1 |  |          |         |          |   m_i   |          |         |          |  |
%  -->O--|          |---...---|          |--- O ---|          |---...---|          |--O m_g
%  u  |  |  |----   |         |  |----   |    |    |  |----   |         |  |----   |  |
%     |  |--| |-----|         |--| |-----|    |    |--| |-----|         |--| |-----|  |
%     |     |----                |----        |       |----                |----      |
%     |      d1_1                   d1_i      |        d1_i+1               d1_g-1    |
%     |                                       |                                       |
%     |      k2_1                   k2_i      |                            k2_g       |
%     |  |-\/\/\/\/-|  |\   /|  |-\/\/\/\/-|  |                    /|  |-\/\/\/\/-|   |
%     |  |          |  |\   /|  |          |  |                    /|  |          |   |
%     ---|          |--|\   /|--|          |---                    /|--|          |----
%        |  |----   |  |\   /|  |  |----   |                       /|  |  |----   |
%        |--| |-----|  |\   /|  |--| |-----|                       /|  |--| |-----|
%           |----                  |----                                  |----
%            d2_1                   d2_i                                   d2_g
%
% REFERENCE:
% [1] V. Mehrmann, T. Stykel. Balanced truncation model reduction for 
%     large-scale systems in descriptor form, in Dimension Reduction of 
%     Large-Scale Systems, P. Benner, V. Mehrmann, and D. Sorensen (eds.),
%     Springer-Verlag, Berlin/Heidelberg, 2005, pp.83--115
%
% ----------------------------------------------------------------------
% T.Stykel, TU Berlin, 9.06.2006

n = 2*g+1;         % state space dimension

% Input parameters are not completely checked
if size(mas,1) < size(mas,2), mas = mas'; end
if size(k1,1)  < size(k1,2),  k1  = k1'; end
if size(k2,1)  < size(k2,2),  k2  = k2'; end
if size(d1,1)  < size(d1,2),  d1  = d1'; end
if size(d2,1)  < size(d2,2),  d2  = d2'; end

% Example from [1]
if length(mas) == 0, mas = 100*ones(g,1); end
if length(k1)  == 0, k1  = 2*ones(g-1,1); end
if length(k2)  == 0, k2  = 2*ones(g,1); k2(1) = 4; k2(g) = 4; end
if length(d1)  == 0, d1  = 5*ones(g-1,1); end
if length(d2)  == 0, d2  = 5*ones(g,1); d2(1) = 10;d2(g) = 10;end

% matrix M
M = spdiags(mas,0,g,g);

% matrix K 
K = spdiags(k1,-1,g,g);
K=K+K'-spdiags([0; k1]+k2+[k1; 0],0,g,g);

% matrix D
D = spdiags(d1,-1,g,g);
D = D+D'-spdiags([0; d1]+d2+[d1; 0],0,g,g);

% matrix G
G = zeros(1,g);
G(1)=1; 
G(end)=-1;

% matrix E
E=sparse(n,n);
E(1:g,1:g)=speye(g);
E(g+1:2*g,g+1:2*g) = M;

% matrix A
A=sparse(n,n);
A(1:g,g+1:2*g)=speye(g);
A(g+1:2*g,1:g)=K;
A(g+1:2*g,g+1:2*g)=D;
A(2*g+1:end,1:g)=G;
A(g+1:2*g,2*g+1:end)=-G';

% matrix B
B=sparse(n,1);
B(g+1,1)=1;

% matrix C
% C=sparse(3,n);
% C(1,1)=1;
% C(2,2)=1;
% C(3,g-1)=1;
C = speye(n);

m=size(B,2); p=size(C,1);
disp('Problem dimensions:');
% nf and ninf are the dimensions of the deflating subspaces of s*E-A
% corresponding to the finite and infinite eigenvalues, n=nf+ninf
disp(['n = ',int2str(n),', nf = ', int2str(n-3), ', ninf = ', int2str(3)]);
disp(['m = ',int2str(m),',  p = ',int2str(p)]);

iM=spdiags(1./mas,0,g,g);  %speye(g)/mas;
G=sparse(1,g); 
G(1,1)=1; G(1,g)=-1;
GG=iM*G'/(G*iM*G');
Pi=speye(g)-GG*G;

% projector Pl
Pl = sparse(n,n);
Pl(1:g,1:g) = Pi;
Pl(1:g,2*g+1:2*g+1) = -Pi*iM*D*GG;
Pl(g+1:2*g,1:g)     = -M*Pi*iM*D*(GG*G);
Pl(g+1:2*g,g+1:2*g) =  M*Pi*iM;
Pl(g+1:2*g,2*g+1)   = -M*Pi*iM*(K+D*Pi*iM*D)*GG;

% projector Pr
Pr = sparse(n,n);
Pr(1:g,1:g)         =  Pi;
Pr(g+1:2*g,1:g)     = -Pi*iM*D*GG*G;
Pr(g+1:2*g,g+1:2*g) =  Pi;
GG = G/(G*iM*G');
Pr(2*g+1,1:g)       =  GG*iM*(K*Pi-D*Pi*iM*D*iM*G'*GG);
Pr(2*g+1,g+1:2*g)   =  GG*iM*D*Pi;