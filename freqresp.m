function [G, omega, sys] = freqresp(sys, varargin)
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

[A,B,C,D,E] = dssdata(sys);
m=sys.m; p=sys.p; n=sys.n;

omegaIndex = cellfun(@isfloat,varargin);
if ~isempty(omegaIndex) && nnz(omegaIndex)
    omega = varargin{omegaIndex};
    varargin{omegaIndex}=[];
else
    [~,omega,sys] = getFreqRange(sys);
end

% iterate over systems in varargin

if (sys.Ts==0) % Convert frequency to either laplace or z variable
    s = 1i* omega;
else
    s = exp(1i* omega*sys.Ts);
end

G=zeros(p,m,length(s));
if not(sys.isBig)
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


function [m,omega,sys] = getFreqRange(sys) %TODO: fix for Ts~=0

% --------- frequency range needs to be chosen ---------
dc = freqrespCell(sys,0);    % G(0)=DCgain
ft = freqrespCell(sys,inf);  % G(inf)=feedthrough

%determine minimum frequency
if any(any(cellfun(@isinf,dc))) || any(any(cellfun(@isnan,dc))) % pole at s=0
    wmin = 0;   %***
    dc = num2cell(ones(size(dc)));
elseif any(any(cellfun(@abs,dc)<1e-14))   % transfer zero at s=0
    wmin = 0;   %***
    dc = num2cell(ones(size(dc)));
else
    wmin=0; t = freqrespCell(sys, 10^wmin);
    while cellfun(@(x,y) norm(x-y)/norm(y),t,dc) > 1e-2
        wmin=wmin-1; t = freqrespCell(sys, 10^wmin);
    end
    while cellfun(@(x,y) norm(x-y)/norm(y),t,dc) < 1e-2
        wmin=wmin+1; t = freqrespCell(sys, 10^wmin);
    end
    wmin=wmin-1;
end

%determine maximum frequency
wmax=0; t = freqrespCell(sys, 10^wmax);
while cellfun(@(x,y,z) norm(x-y)/norm(z),t,ft,dc) > 1e-6
    wmax=wmax+1; t = freqrespCell(sys, 10^wmax);
end
while cellfun(@(x,y,z) norm(x-y)/norm(z),t,ft,dc) < 1e-6
    wmax=wmax-1; t = freqrespCell(sys, 10^wmax);
end
wmax=wmax+1;

delta = (wmax-wmin)/19; % initial resolution (insert odd number only!)
omega = 10.^(wmin:delta:wmax);
m = freqrespCell(sys, omega);

while(1)
    % increase plot density until vertical change per step is below 1%
    for k=1:length(omega)-1
        if cellfun(@(x) abs(abs(x(k)) - abs(x(k+1)))/(abs(x(k)) + abs(x(k+1))),m) > 0.01
            break
        end
    end
    if k==length(omega)-1
        break
    end
    % do not refine above 2000 points
    if length(omega)>1000
        break
    end
    delta = delta/2;
    omega = 10.^(wmin:delta:wmax);
    
    % calculate new values of frequency response
    temp=freqrespCell(sys, omega(2:2:length(omega)));
    
    % update array of results (insert new values)
    m=cellfun(@(x,y) [reshape([x(1:length(x)-1);y],1,2*length(x)-2),x(end)],m,temp,'UniformOutput',false);
end

% determine magnitude and phase from complex frequency response
mag = cellfun(@abs, m, 'UniformOutput',false);

% determine H_inf-norm (maximum magnitude)
[a,b]=cellfun(@max, mag);
[a,c]=max(a);
[H_inf_norm,d]=max(a);
H_inf_peakfreq=omega(b(c(d),d));
% Issue: the values get assigned but this does not affect the original
% sss object as it is called by value not by reference as @sss is not a
% handle class.
sys.H_inf_norm = max([sys.H_inf_norm, H_inf_norm]);
sys.H_inf_peakfreq = max([sys.H_inf_peakfreq, H_inf_peakfreq]);

end

function [m, omega] = freqrespCell(varargin)
[m, omega] = freqresp(varargin{:});
warning('off', 'MATLAB:mat2cell:TrailingUnityVectorArgRemoved');
m = mat2cell(m,ones(size(m,1),1),ones(size(m,2),1),size(m,3));
m = cellfun(@(x) x(:,:),m, 'UniformOutput', false);
end