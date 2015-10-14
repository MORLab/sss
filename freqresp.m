function G = freqresp(sys, s, varargin)
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
%              Lisa Jeschek
% ------------------------------------------------------------------

[A,B,C,D,E] = ABCDE(sys);
m=sys.m; p=sys.p; n=sys.n;
if nargin >= 4
    % I/O-pair selected
    if ~isempty(varargin{1}) && ~isnan(varargin{1})
        B=B(:,varargin{1});
        D=D(:,varargin{1});
        m=size(B,2);
    end
    if ~isempty(varargin{2}) && ~isnan(varargin{2})
        C=C(varargin{2},:);
        D=D(varargin{2},:);
        p=size(C,1);
    end
end

G=zeros(p,m,length(s));
if size(A,1) < 5000
    for i=1:length(s)
        G(:,:,i) = freqresp_local(A,B,C,D,E, s(i),n);
    end
else
    parfor i=1:length(s)
        G(:,:,i) = freqresp_local(A,B,C,D,E, s(i),n);
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
