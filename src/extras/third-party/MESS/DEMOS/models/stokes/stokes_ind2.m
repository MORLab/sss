function [ E, A, B, C, nf ] = stokes_ind2(m, q, nx, ny)
%
% Semidiscretized 2D Stokes equation:
%
%       E*x'(t) = A*x(t) + B*u(t),
%          y(t) = C*x(t),
%
% with
%      E = [ I  0 ], A = [ A11 A12 ], B = [B1], C =[ C1, C2 ],
%          [ 0  0 ]      [ A21  0  ]      [B2]
% where
%   A11 is the discrete Laplace operator,
%  -A12 is the discrete gradient operator,
%  -A21 is the discrete divergence operator (here A21=A12^T),
%   ( A12 and A21' have full column rank => s*E-A is of index 2 ).
%
% INPUT:
%       m  number of inputs
%       q  number of outputs
%      nx  number of grid points in the x-direction ( nx>1 )
%      ny  the number of grid points in the y-direction ( ny>1 )
%
% OUTPUT:
%      E   real n-by-n sparse matrix with
%           n = (nx-1)*ny+(ny-1)*nx+nx*ny-1
%      A   real n-by-n sparse matrix
%      B   real n-by-m matrix, m is the number of inpits
%      C   real q-by-n sparse matrix, q is the number of outputs
%     nf   dimension of the deflating subspaces corresponding to
%            the finite eigenvalues
%
% REFERENCE:
%  M.Schmidt. Systematic discretization of input/output maps and
%  other contributions to the control of distributed parameter systems.
%  Ph.D. thesis, TU Berlin, 2007.
%  http://dx.doi.org/10.14279/depositonce-1600
%  (Section 3.3.4, pp. 34-38, Section 3.7.1, pp. 62-63)
% ----------------------------------------------------------------------
% M.Schmidt, TU Berlin, 31.05.2007
% T.Stykel, TU Berlin, 24.11.2007


global OmegaIn OmegaOut mu_choice nu_choice

%----------------------------------------------------------
% parameters

% geometry: rectangular domain
Length1 = 1;
Length2 = 1;

% geometry: rectangular control and observability subdomains
OmegaIn     = [0.1,0.9,0.1,0.3];
OmegaOut    = [0.4,0.6,0.4,0.9];

% choose input and output basis functions
mu_choice   = 1; % see bfnc_mu.m for meaning of 1
nu_choice   = 1;

%-----------------------------------------------------------
% initialization

% some often used variables...
h1      = Length1/nx;
h2      = Length2/ny;
x1      = linspace(h1/2,Length1-h1/2,nx);
x1e1    = linspace(h1,Length1-h1,nx-1);
x2      = linspace(h2/2,Length2-h2/2,ny);
x2e2    = linspace(h2,Length2-h2,ny-1);

nv1 = (nx-1)*ny;
nv2 = (ny-1)*nx;
nv  = nv1+nv2;
np  = nx*ny;
nf = (nx-1)*(ny-1);

% preallocate matrices
L1  = sparse(nv1,nv1);
G1  = sparse(nv1,np);
B1  = sparse(nv1,m);
C1  = sparse(q,nv1);
D1  = sparse(np,nv1);

L2  = sparse(nv2,nv2);
G2  = sparse(nv2,np);
B2  = sparse(nv2,m);
C2  = sparse(q,nv2);
D2  = sparse(np,nv2);

fprintf(1,'degrees of freedom: \n');
fprintf(1,'------------------------------\n');
fprintf(1,' total         :%6.0f\n',nv1+nv2+np);
fprintf(1,' velocity      :%6.0f\n',nv1+nv2);
fprintf(1,' pressure      :%6.0f\n',np);
fprintf(1,' n_finite      :%6.0f\n',nf);
fprintf(1,'------------------------------\n');

info.domain     = [0,Length1,0,Length2];
info.FVMdisc    = [nx,ny] ;
info.OmegaIn    = OmegaIn;
info.OmegaOut   = OmegaOut;
info.Ndof.total = nv1+nv2+np;
info.Ndof.v1    = nv1;
info.Ndof.v2    = nv2;
info.Ndof.pres  = np;
info.Ndof.nf    = nf;

%------------------------------------------------------------
% build FVM matrices
disp('Generating FVM matrices...');

% matrices divided by (h1*h2) in order to have mass matrix=identity
% Matrix L1
disp(' -> Laplacians...');
auxM = sparse(nx-1,nx-1);
auxM = auxM - 2* speye(nx-1) + sparse(diag(ones(nx-2,1),1));
auxM = auxM + sparse(diag(ones(nx-2,1),-1));
L1   = L1 + 1/(h1*h1) * kron(speye(ny),auxM);
clear auxM;

% Matrix L2
auxM = sparse(ny-1,ny-1);
auxM = auxM - 2* speye(ny-1) + sparse(diag(ones(ny-2,1),1));
auxM = auxM + sparse(diag(ones(ny-2,1),-1));
L2   = L2 + 1/(h2*h2) * kron(auxM,eye(nx));
clear auxM;

% Matrices G1 and D1
disp(' -> gradient and divergence operator...');
auxM = sparse(nx-1,nx);
nVec = sparse(nx-1,1);
auxM = auxM + [nVec,speye(nx-1)] - [speye(nx-1),nVec];
G1   = G1 + 1/h1 * kron(speye(ny),auxM);
D1   = -G1';
clear auxM;

% Matrices G2 and D2
auxM = sparse(ny-1,ny);
nVec = sparse(ny-1,1);
auxM = auxM + [nVec,speye(ny-1)] - [speye(ny-1),nVec];
G2   = G2 + 1/h2 * kron(auxM,eye(nx));
D2   = -G2';
clear auxM;

% B1 and v1 velocity nodes
disp(' -> B1, C1 and v1 velocity nodes...');
B1=sparse(nv1,m);
C1=sparse(q,nv1);
for i2=1:ny,
    xx1 = x1e1;
    xx2 = x2(i2)*ones(1,nx-1);
    v1x1((i2-1)*(nx-1)+1:i2*(nx-1),1)=xx1;
    v1x2((i2-1)*(nx-1)+1:i2*(nx-1),1)=xx2;
    for k=1:m,
        B1((i2-1)*(nx-1)+1:i2*(nx-1),k)=vecBmu1(xx1,xx2,k)';
    end;
    for l=1:q,
        C1(l,(i2-1)*(nx-1)+1:i2*(nx-1))=h1*h2*vecCstarnu1(xx1,xx2,l);
    end;
end;

% B2 and v2 velocity nodes
disp(' -> B2, C2 and v2 velocity nodes...');
B2=sparse(nv2,m);
C2=sparse(q,nv2);
for i2=1:ny-1,
    xx1 = x1;
    xx2 = x2e2(i2)*ones(1,nx);
    v2x1(1,(i2-1)*nx+1:i2*nx)=xx1;
    v2x2(1,(i2-1)*nx+1:i2*nx)=xx2;
    for k=1:m,
        B2((i2-1)*nx+1:i2*nx,k)=vecBmu2(xx1,xx2,k)';
    end;
    for l=1:q,
        C2(l,(i2-1)*nx+1:i2*nx)=h1*h2*vecCstarnu2(xx1,xx2,l);
    end;
end;

%-----------------------------------------------------
% extract DAE

disp('Setting up system matrices ...');
E    = sparse(nv1+nv2+np-1,nv1+nv2+np-1);
A    = sparse(nv1+nv2+np-1,nv1+nv2+np-1);
B    = sparse(nv1+nv2+np-1,m);
C    = sparse(q,nv1+nv2+np-1);
v1xx = [v1x1;v1x2];
v2xx = [v2x1;v2x2];
info = info;

%E(1:nv1+nv2,1:nv1+nv2)         = speye(nv1+nv2);
for j=1:nv1+nv2
    E(j,j)=1;
end

A = [L1 sparse(nv1,nv2) -G1(:,1:np-1);
     sparse(nv2,nv1) L2 -G2(:,1:np-1);
     D1(1:np-1,:) D2(1:np-1,:) sparse(np-1,np-1)];

B(1:nv1+nv2,:)  = [B1;B2];
C(:,1:nv1+nv2)  = [C1,C2];

% % projectors Pl and Pr
% A12=[-G1(:,1:np-1); -G2(:,1:np-1)];
% AA=A12'*A12;
% Pi=speye(nv)-A12*(AA\A12');
% Pl = [Pi -Pi*(A(1:nv,1:nv)*((AA\(A12'))')); sparse(np-1,nv+np-1)];
% Pr = [Pi sparse(nv,np-1); -AA\(A12'*(A(1:nv,1:nv)*Pi)) sparse(np-1,np-1)];

clear L1 L2 G1 G2 D1 D2;

%-------------------------------------------------------
function [out,info]=Bmu(x1,x2,bfi,t);
%  Bmu(x1,x2,bfi)
%
% Version 1.2, last change: 25-07-06
% Michael Schmidt, TU Berlin, Inst. f. Mathematik,
% mschmidt@math.tu-berlin.de
%
% out = Bmu(x1,x2,bfi,t)

global controlChoice OmegaIn;

% specify default input function for the case nargin==3
defaultInputName    ='u';

% specify support rectangle
if isempty('OmegaIn'), error('Bmu.m requires global "OmegaIn"'); end;
         dc1 = (OmegaIn(2)-OmegaIn(1))/2;
         xc1 = OmegaIn(1)+dc1;
         dc2 = (OmegaIn(4)-OmegaIn(3))/2;
         xc2 = OmegaIn(3)+dc2;
info.Omega_c = OmegaIn;

if size(x1,1)>1, x1=x1'; end;
if size(x2,1)>1, x2=x2'; end;

% preallocate output
out=zeros(size(x1));

% find points in \Omega_c
PInd=intersect(find(x1>xc1-dc1 & x1<xc1+dc1) , find(x2>xc2-dc2 & x2<xc2+dc2));

% theta(x1)
theta   =(x1(PInd)-xc1)/(2*dc1)+0.5;  %x1 -->theta \in [0,1]
x2tilde =(x2(PInd)-xc2)/(2*dc2)+0.5;  %x2 -->x2tilde \in [0,1]

% calculate Bmu
switch nargin
    case 3, % Bmu as initial value or observation weight with mu-basis funtion
        bfi=bfi(1);
        out(PInd)=bfnc_mu(theta,bfi).*w_c(x2tilde);
        % specify support of Bmu
        [forget,mu_info]=bfnc_mu(0,bfi);
        info.suppBmu =[ (mu_info.suppstart-0.5)*2*dc1+xc1 (mu_info.suppend-0.5)*2*dc1+xc1  xc2-dc2 xc2+dc2];
    case 4, %control right hand sight for control function u(theta,t)
        t=t(1);
        % specify control input function
        if exist('controlChoice'),
            inputName = controlChoice;
            if isempty(controlChoice);
                warning(['default input used: "',defaultInputName,'.m"']);
                inputName=defaultInputName;end;
        else,
            inputName   =defaultInputName;
        end;
        % u(theta,t) has dimensions length(theta) x length(t)
        out(PInd)=feval(inputName,theta,t)'.*w_c(x2tilde);
    otherwise
        error('not right choice of input arguments!');
end;

%------------------
% weighting function w_c:[0,1] --> R
function out=w_c(x2);
global Wn
Wn=[1];
out=zeros(size(x2));
for n=1:length(Wn),
    out=out+sin(n*pi*x2);
end;

%---------------------------------------------------------
function [out,info]=vecBmu1(x1,x2,bfi);

bfi=bfi(1); %necessary since size(bfi)=size(x1);
out=zeros(size(x1));

switch nargin,
    case 3,
        if mod(bfi,2)==1,
            [out,info]=Bmu(x1,x2,(bfi+1)/2);
        end;
    otherwise,
        error('false number of input arguments');
end;
out=100*out;

%---------------------------------------------------------
function [out,info]=vecBmu2(x1,x2,bfi);

bfi=bfi(1);%necessary since size(bfi)=size(x1);
out=zeros(size(x1));

switch nargin,
    case 3,
        if mod(bfi,2)==0,
            [out,info]=Bmu(x1,x2,bfi/2);
        end;
    otherwise,
        error('false number of input arguments');
end;
out=100*out;

%--------------------------------------------------------
function [out,mu_info]=bfnc_mu(x,bfi,lcomb)
% [out,mu_info]=inpbasis(x,bfi)
% creates continuous basis functions \mu_i of \tilde U
% output is always row vector
%
% depends on global variable mu_basis with mu_basis = ...
% 1 : linear hat functions in L^2(0,1):   hierarchical
% 2 : linear hat functions in H_0^1(0,1): hierarchical
% 3 : linear hat functions in L^2(0,1):   nodal
% 4 : linear hat functions in H_0^1(0,1): nodal
% 5 : sinus/cosinus in L^2(0,1)
% 6 : sinus
%
% mu_info.suppstart  - begin of support
% mu_info.suppend    - end of support
% mu_info.node       - position of node
% mu_info.bfi_lvl    - level of specific basis function bfi
% mu_info.bfi_lvlind - index on this level
% mu_info.h1l        - grid fineness on this level

global mu_choice pmax

%disp('bnfc_mu');

if nargin<2, error('not enough input arguments'); end;
if isempty(mu_choice), error('mu_choice must be globally defined');end;


if nargin==2,
    switch mu_choice,
        case 2, %linear hat functions in H_0^1: hierarchical
            % find level and levelindex
            lvl     =length(dec2bin(bfi))-1;
            lvlind  =bfi-2^lvl;
            % scale and shift reference function
            out=PSI2(2^lvl*x-lvlind*ones(size(x)));
            % give information on basis function
            mu_info.suppstart  = lvlind/2^lvl;
            mu_info.node       = (lvlind+0.5)/2^lvl;
            mu_info.suppend    = (lvlind+1)/2^lvl;
            %mu_info.bflevels   = 0;
            mu_info.bfi_lvl    = lvl;
            mu_info.bfi_lvlind = lvlind;
            mu_info.h          = 1/2^(lvl+1);
        case 1, %linear hat functions in L^2: hierarchical
            if bfi<=2, %two additional basis functions w.r.t. case 1
                if bfi==1,
                    out=ones(size(x))-x;
                    mu_info.node=0;
                else %bfi==2
                    out=x;
                    mu_info.node=1;
                end;
                % give information on basis function
                mu_info.suppstart  = 0;
                mu_info.suppend    = 1;
                mu_info.node       = [];
                mu_info.bfi_lvl    = -1;
                mu_info.bfi_lvlind = bfi;
                mu_info.h          = 1;
            else, % as in case 1 with shifted bfi
                bfi=bfi-2;
                % find level and levelindex
                lvl     =length(dec2bin(bfi))-1;
                lvlind  =bfi-2^lvl;
                % scale and shift reference function
                out=PSI2(2^lvl*x-lvlind*ones(size(x)));
                % give information on basis function
                mu_info.suppstart  = lvlind/2^lvl;
                mu_info.node       = (lvlind+0.5)/2^lvl;
                mu_info.suppend    = (lvlind+1)/2^lvl;
                %mu_info.bflevels   = 0;
                mu_info.bfi_lvl    = lvl;
                mu_info.bfi_lvlind = lvlind;
                mu_info.h          = 1/2^(lvl+1);
            end;
        case 3,
            if isempty(pmax), error('global variable pmax not defined...');end;
            h=1/(pmax-1);
            out=PSI2(x/(2*h)-(bfi-2)/2);
            mu_info.suppstart  = max(0,(bfi-2)*h);
            mu_info.node       = (bfi-1)*h;
            mu_info.suppend    = min(1,bfi*h);
            mu_info.h          =h;
        case 4,
            error('not yet implemented');
        case 5,
            if mod(bfi,2)==1, %impair bfi --> cosinus (cos(0) is first basis function!)
                if floor(bfi/2)==0,
                    out=ones(size(x));
                else,
                    out=sqrt(2)*cos(2*floor(bfi/2)*pi*x);
                end;
            else, %pair bfi --> sinus
                out=sqrt(2)*sin(2*bfi/2*pi*x);
            end
            mu_info.suppstart  = 0;
            mu_info.node       = [];
            mu_info.suppend    = 1;

        case 6, %FOR COMPARISON WITH EXACT KERNEL (HOM. DIRICHLET ON [0,1]^2)
            out=sqrt(2)*sin(bfi*pi*x);
            mu_info.suppstart  = 0;
            mu_info.node       = [];
            mu_info.suppend    = 1;

        otherwise
            error('Basis or method choice invalid');
    end;
else, % LINEAR COMBINATION lcomb
    out = zeros(size(x));
    mu_info.suppstart   = 1;
    mu_info.node        = [];
    mu_info.suppend     = 0;
    for n=1:length(lcomb),
        [out_temp,info]     = bfnc_mu(x,n);
        out                 = out+lcomb(n)*out_temp;
        mu_info.suppstart   = min(mu_info.suppstart,info.suppstart);
        mu_info.suppend     = max(mu_info.suppend,info.suppend);
    end;
end;

if size(out,1)>1, out=out'; end;

%--------------------------------------------------------
function [out,nu_info]=bfnc_nu(x,bfi,lcomb)

global nu_choice

% [out,nu_info]=inpbasis(x,bfi)
% creates continuous basis functions \nu_i of \tilde Y
%
% depends on global variable nu_basis with nu_basis = ...
% 1 : linear hat functions in L^2(0,1):   hierarchical
% 2 : linear hat functions in H_0^1(0,1): hierarchical
% 3 : linear hat functions in H_0^1(0,1): nodal
% 4 : linear hat functions in L^2(0,1):   nodal
% 5 : sinus/cosinus in L^2(0,1)
%
% nu_info.suppstart  - begin of support
% nu_info.suppend    - end of support
% nu_info.node       - position of node
% nu_info.bfi_lvl    - level of specific basis function bfi
% nu_info.bfi_lvlind - index on this level
% nu_info.h1l        - grid fineness on this level

if nargin<2, error('not enough input arguments'); end;
if isempty(nu_choice), error('nu_choice must be globally defined');end;

if nargin==2,
switch nu_choice,
    case 2, %linear hat functions in H_0^1: hierarchical
        % find level and levelindex
        lvl     =length(dec2bin(bfi))-1;
        lvlind  =bfi-2^lvl;
        % scale and shift reference function
        out=PSI2(2^lvl*x-lvlind*ones(size(x)));
        % give information on basis function
        nu_info.suppstart  = lvlind/2^lvl;
        nu_info.node       = (lvlind+0.5)/2^lvl;
        nu_info.suppend    = (lvlind+1)/2^lvl;
        %nu_info.bflevels   = 0;
        nu_info.bfi_lvl    = lvl;
        nu_info.bfi_lvlind = lvlind;
        nu_info.h          = 1/2^(lvl+1);
    case 1, %linear hat functions in L^2: hierarchical
        if bfi<=2, %two additional basis functions w.r.t. case 1
            if bfi==1,
                out=ones(size(x))-x;
                mu_info.node=0;
            else %bfi==2
                out=x;
                nu_info.node=1;
            end;
            % give information on basis function
            nu_info.suppstart  = 0;
            nu_info.suppend    = 1;
            nu_info.bfi_lvl    = -1;
            nu_info.bfi_lvlind = bfi;
            nu_info.h          = 1;
        else, % as in case 1 with shifted bfi
            bfi=bfi-2;
            % find level and levelindex
            lvl     =length(dec2bin(bfi))-1;
            lvlind  =bfi-2^lvl;
            % scale and shift reference function
            out=PSI2(2^lvl*x-lvlind*ones(size(x)));
            % give information on basis function
            nu_info.suppstart  = lvlind/2^lvl;
            nu_info.node       = (lvlind+0.5)/2^lvl;
            nu_info.suppend    = (lvlind+1)/2^lvl;
            %nu_info.bflevels   = 0;
            nu_info.bfi_lvl    = lvl;
            nu_info.bfi_lvlind = lvlind;
            nu_info.h          = 1/2^(lvl+1);
        end;
    case 3,
        error('not yet implemented');
    case 4,
        error('not yet implemented');
    case 5,
        if mod(bfi,2)==1, %impair bfi --> cosinus (cos(0) is first basis function!)
            out=0.5*cos(2*floor(bfi/2)*pi*x);
        else, %pair bfi --> sinus
            out=0.5*sin(2*bfi/2*pi*x);
        end
        nu_info.suppstart  = 0;
        nu_info.suppend    = 1;
    otherwise
        error('Basis or method choice invalid');
end;
else,
        out=zeros(size(x));
    nu_info.suppstart=1;
    nu_info.suppend  =0;
    for n=1:length(lcomb),
        [out_temp,info]=bfnc_nu(x,n);
        out=out+lcomb(n)*out_temp;
        nu_info.suppstart=min(nu_info.suppstart,info.suppstart);
        nu_info.suppend  =max(nu_info.suppend,info.suppend);
    end;
end;

%---------------------

% linear hat function with support [0,1]
function out=PSI2(t);
ind_p1=find(t>=0&t<0.5);
ind_m1=find(t>=0.5&t<1);
out=zeros(size(t));
out(ind_p1)=t(ind_p1)*2;
out(ind_m1)=2-t(ind_m1)*2;

%-------------------------------------------------------
function [out,info]=Cstarnu(x1,x2,bfi);
%
% Version 1.2, last change: 25-07-06
% Michael Schmidt, TU Berlin, Inst. f. Mathematik,
% mschmidt@math.tu-berlin.de

% initialization
global OmegaOut
bfi     =bfi(1);
xi_of   ='x2';
out     =zeros(size(x1));

% specify support rectangle
if isempty('OmegaOut'), error('Cstarnu.m requires global "OmegaOut"'); end;
dm1=(OmegaOut(2)-OmegaOut(1))/2;
xm1=OmegaOut(1)+dm1;
dm2=(OmegaOut(4)-OmegaOut(3))/2;
xm2=OmegaOut(3)+dm2;
info.Omega_m =OmegaOut;

switch xi_of
    case 'x1',
        % find points in \Omega_m and assign value there
        fInd=intersect(find(x1>xm1-dm1 & x1<xm1+dm1) , find(x2>xm2-dm2 & x2<xm2+dm2));

        xi          =(x1(fInd)-xm1)/(2*dm1)+0.5;
        out(fInd)   =bfnc_nu(xi,bfi)/(2*dm2*2*dm1);

        % specify support of Cstarnu
        [forget,nu_info]   =bfnc_nu(0,bfi);
        info.suppCstarnu   =[(nu_info.suppstart-0.5)*2*dm1+xm1...
                            (nu_info.suppend-0.5)*2*dm1+xm1...
                            xm2-dm2 xm2+dm2];
    case 'x2',
        % find points in \Omega_m and assign value there
        fInd=intersect(find(x1>xm1-dm1 & x1<xm1+dm1) , find(x2>xm2-dm2 & x2<xm2+dm2));

        xi          =(x2(fInd)-xm2)/(2*dm2)+0.5;
        out(fInd)   =bfnc_nu(xi,bfi)/(2*dm2*2*dm1);

        % specify support of Cstarnu
        [forget,nu_info]   =bfnc_nu(0,bfi);
        info.suppCstarnu   =[xm1-dm1 xm1+dm1...
                            (nu_info.suppstart-0.5)*2*dm2+xm2...
                            (nu_info.suppend-0.5)*2*dm2+xm2];
end;

%-----------------------------------------------------
function [out,info]=vecCstarnu1(x1,x2,bfi);

bfi=bfi(1);
out=zeros(size(x1));
switch nargin,
    case 3,
        if mod(bfi,2)==1,
            [out,info]=Cstarnu(x1,x2,(bfi+1)/2);
        end;
    otherwise,
        error('false number of input arguments');
end;

%-----------------------------------------------------
function [out,info]=vecCstarnu2(x1,x2,bfi);

bfi=bfi(1);
out=zeros(size(x1));
switch nargin,
    case 3,
        if mod(bfi,2)==0,
            [out,info]=Cstarnu(x1,x2,bfi/2);
           % disp(['calling Cstarnu ',numstr(bfi/2)]);
        end;
    otherwise,
        error('false number of input arguments');
end;