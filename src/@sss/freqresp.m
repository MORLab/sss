function [G, omega, sys] = freqresp(varargin)
% freqresp - Evaluates complex transfer function of LTI systems
% 
% Description:
%       Evaluates complex transfer function of LTI systems. If the vector
%       of complex frequencies is not passed, then a range of imaginary
%       frequencies is automatically selected.
%
% Syntax:
%       G = freqresp(sys, s)
%       G = freqresp(sys, s, opts)
%       [G, omega] = freqresp(sys)
%
% Inputs:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       -s: vector of complex frequencies
%       *Optional Input Arguments:*
%       -opts: Computation options, see <a href="matlab:help freqresp">freqresp</a>
%       
% Outputs:      
%       -G: vector of complex frequency response values
%
% Examples:
%
% See Also:
%
% References:
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Chair of Thermofluid Dynamics, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Stefan Jaensch, Heiko Panzer, Sylvia Cremer, Rudy Eid
%               Lisa Jeschek, Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  07 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% input parsing
sys= varargin{1};

omegaIndex = cellfun(@isfloat,varargin);
if ~isempty(omegaIndex) && nnz(omegaIndex)
    %frequency vector was specified
    omega = varargin{omegaIndex};
    varargin(omegaIndex)=[];
else
    %frequency vector will be identified
    [G,omega,sys] = getFreqRange(sys);
    return
end

%%  Compute the value of the transfer function at selected freq.
if (sys.Ts==0) % Convert frequency to either laplace or z variable
    s = 1i* omega;
else
    s = exp(1i* omega*sys.Ts);
end

for iSys=1:length(varargin)
    sys= varargin{iSys};
    m=sys.m; p=sys.p; n=sys.n;
    [A,B,C,D,E] = dssdata(sys);
    G=zeros(p,m,length(s));
    for iShift=1:length(s)
        G(:,:,iShift) = freqresp_local(A,B,C,D,E, s(iShift),n);
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

function [m,omega,sys] = getFreqRange(sys)

% frequency range needs to be chosen
dc = freqresp(sys,0);    % G(0)=DCgain
ft = freqresp(sys,inf);  % G(inf)=feedthrough

%determine minimum frequency
if any(any(isinf(dc))) || any(any(isnan(dc))) % pole at s=0
    wmin = 0;  dc = ones(size(dc));
elseif any(any(abs(dc<1e-14)))   % transfer zero at s=0
    wmin = 0;  dc = ones(size(dc));
else
    wmin=0; t = freqresp(sys, -1i*10^wmin);
    while norm(t-dc)/norm(dc) > 1e-2, 
        wmin=wmin-1; t = freqresp(sys, 10^wmin);
    end
    while norm(t-dc)/norm(dc) < 1e-2
        wmin=wmin+1; t = freqresp(sys, 10^wmin);
    end
    wmin=wmin-1;
end

%determine maximum frequency
wmax=0; t = freqresp(sys, -1i*10^wmax);
while norm(t-ft)/norm(dc) > 1e-6
    wmax=wmax+1; t = freqresp(sys, 10^wmax);
end
while norm(t-ft)/norm(dc) < 1e-6
    wmax=wmax-1; t = freqresp(sys, 10^wmax);
end
wmax=wmax+1;

delta = (wmax-wmin)/19; % initial resolution (insert odd number only!)
omega = 10.^(wmin:delta:wmax);
m = freqresp(sys, omega);

while(1)
    % increase plot density until vertical change per step is below 1%
    if ~any(any(abs(diff(m,1,3))./abs(m(:,:,1:end-1)) > 0.01))
        % relative accuracy achieved
        break
    end
    if length(omega)>1000 % do not refine above 1000 points
        break
    end
    delta = delta/2;
    omega = 10.^(wmin:delta:wmax);
    
    % calculate new values of frequency response
    temp=freqresp(sys, omega(2:2:length(omega)));
    
    % update array of results (insert new values)
    mOld = m; m = zeros(size(m,1),size(m,2),length(omega));
    m(:,:,2:2:length(omega)) = temp; m(:,:,1:2:end) = mOld;
end

end
