function sysd = diag(sys)
% DIAG - Transforms an LTI system to (block)diagonal representation
%
% Syntax:
%       sysd = DIAG(sys)
%
%
% Description:
%       sysd = DIAG(sys) transforms the sparse state-space model sys to a 
%       (block)diagonal representation. TODO
%       
%       During the diagonalization, the C-vector is normalized to contain
%       ones.
%
%
% Input Arguments:
%       -sys: an sss-object containing the LTI system
%
%
% Output Arguments:
%       -sysd: system in diagonal form
%
%
% Examples:
%       To compute the diagonal state-space realization of the benchmark
%       "build" (SISO) use
%
%>      load build.mat
%>      sys = sss(A,B,C)
%>      sysd = diag(sys)
%>      [Ad,Bd,Cd] = ssdata(sysd);
%
%>      figure;
%>      subplot(1,2,1); spy(sys.A); title('Sparsity pattern of A');
%>      subplot(1,2,2); spy(Ad); title('Sparsity pattern of Ad');
%
%       DIAG also supports SIMO, MISO and MIMO systems:
% 
%>      load PEEC_MTLn1600.mat
%>      sys = sss(A,B,C,zeros(14,14),E)
%>      sysd = diag(sys)
%>      [Ad,Bd,Cd,Dd] = ssdata(sysd);
%
%
% See also: 
%		eig, residue
%
%
% ------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
% ------------------------------------------------------------------
% Authors:      Heiko K.F. Panzer, Maria Cruz Varona
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  29 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%perform eigen-decomposition of system
if sys.isDescriptor
    [T,A] = eig(full(sys.A),full(sys.E));
    if max(max(A))==inf %isDae == 1
        error('System contains algebraic states (DAE).')
    end
else
    [T,A] = eig(full(sys.A));
end

% transform system to diagonal form
B = (sys.E*T)\sys.B;
C = sys.C*T;

% save poles and residues
sys.poles = transpose(diag(A));
if sys.isSimo || sys.isMiso || sys.isMimo
    %simo, miso or mimo
    r = cell(sys.p, sys.m);
    for j=1:sys.m
        for i=1:sys.p
            r{i,j} = C(i,:) .* conj(B(:,j)');
        end
    end
else
    %siso
    r = {C .* conj(B')};
end
sys.residues = r;

% find real system representation
i = 0;
while i<length(B)
    i=i+1;
    if i<length(B)
        if abs(real(A(i,i)+A(i+1,i+1))) >= abs(imag(A(i,i)+A(i+1,i+1)))*10^3 && ...
           abs(imag(A(i,i)-A(i+1,i+1)))/abs(real(A(i,i)+A(i+1,i+1))) >= 10^(-3)
            delta = real(A(i,i) + A(i+1,i+1))/2;
            omega = imag(A(i,i) - A(i+1,i+1))/2;
            A(i:i+1,i:i+1) = [delta, omega; -omega, delta];
            
            % residues are shifted to input vector
            pc = real(B(i)*C(i)+B(i+1)*C(i+1))/2;
            pd = imag(B(i)*C(i)-B(i+1)*C(i+1))/2;
            B(i) = pc+pd;
            B(i+1) = pc-pd;
            i=i+1;
            continue
        end
    end
    B(i) = B(i) * C(i);
end

% remove remaining imaginary components (resulting from numerical noise)
B = real(B);

sysd = sss(A, B, ones(size(C)), sys.D);

if inputname(1)
    assignin('caller', inputname(1), sys);
end
