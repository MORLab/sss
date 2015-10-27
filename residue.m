function [r,p,d] = residue(sys)
% Computes residues, poles and feedthrough of an LTI system
% ------------------------------------------------------------------
% [r,p,d] = residue(sys)
% Inputs:       * sys: an sss-object containing the LTI system
% Outputs:      * residues r, poles p and feedthrough d, such that
%                   G(s) = r_i/(p_i+s) + d
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer
% Last Change:  19 Out 2015
% ------------------------------------------------------------------

% are poles and residues already available?
if ~isempty(sys.poles) && ~isempty(sys.residues)
    p=sys.poles;
    r=sys.residues;
    d=sys.D;
    return
end

%perform eigen-decomposition of system
try
    [T,J] = eig(sys);
catch err
    error('Computation of the eigenvalues and eigenvectors failed with message:%s',err.message);
end

% transform system to diagonal form
p=diag(J).';
if issparse(T)
    rcondNumber = 1/condest(T);
else
    rcondNumber=rcond(T);
end
if rcondNumber<eps
    warning('Matrix of eigenvectors is close to singular or badly scaled. Results may be inaccurate. RCOND =',num2str(rcondNumber));
end
B=(sys.E*T)\sys.B;
C=sys.C*T;
d=sys.D;

% calculate residues
if sys.isMimo
    %mimo
    r = cell(sys.n, 1);
    for i=1:sys.n
        r{i} = C(:,i)*B(i,:);
    end
else
    %siso
    r = {C .*B.'};
end

% store results to caller workspace
sys.residues=r;
sys.poles=p;
if inputname(1)
    assignin('caller', inputname(1), sys);
end
return
