function G = freqresp(sys, s)
% Evaluates complex transfer function of LTI systems
% ------------------------------------------------------------------
% G = freqresp(sys, s, varargin)
% Inputs:       * sys: an sss-object containing the LTI system
%               * s: vector of complex frequencies
% Outputs:      * G: vector of complex frequency response values
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:     Stefan Jaensch, Heiko Panzer (heiko@mytum.de), Sylvia Cremer, Rudy Eid
% ------------------------------------------------------------------

[A,B,C,D,E] = ABCDE(sys);
m=sys.m; p=sys.p;

G=zeros(p,m,length(s));
if size(A,1) < 5000
    for i=1:length(s)
        G(:,:,i) = freqresp_local(A,B,C,D,E, s(i),sys.n);
    end
else
    parfor i=1:length(s)
        G(:,:,i) = freqresp_local(A,B,C,D,E, s(i),sys.n);
    end
end
end

function G = freqresp_local(A,B,C,D,E, s,n)
% calculate value of transfer function for given vector of frequencies

if isinf(s) || n==0
    G = D;
    return
end

[L,U,k,l,S]=lu(A - s*E, 'vector');
warning('OFF', 'MATLAB:nearlySingularMatrix')
b=S\B; b=b(k,:);
x=L\b;
x(l,:)=U\x;
warning('ON', 'MATLAB:nearlySingularMatrix')
G = -C*x + D;
end
