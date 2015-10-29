function  [varargout] = bode(varargin)
% Plots the bode diagram of an LTI system
% ------------------------------------------------------------------
% [mag, phase, omega] = bode(sys, omega, in, out, options)
% Inputs:       * sys: an sss-object containing the LTI system
%    [optional] * vector of imaginary frequencies to plot at
%               * plot options. see <a href="matlab:help plot">PLOT</a>
% Outputs:      * vector of complex frequency response values
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Stefan Jaensch,
%               Sylvia Cremer, Rudy Eid, Alessandro Castagnotto,
%               Lisa Jeschek
% Last Change:  14 Oct 2015
% ------------------------------------------------------------------


% --------------- EVALUATE OPTIONS ---------------
omegaIndex = cellfun(@isfloat,varargin);
if ~isempty(omegaIndex) && nnz(omegaIndex)
    omega = varargin{omegaIndex};
    varargin{omegaIndex}=[];
end

for i = 1:length(varargin)
    if isa(varargin{i},'sss')
        nSys=i;
    end
end

% --------------- CALCULATE FREQUENCY RESPONSE ---------------
for i = 1:length(varargin)
    if isa(varargin{i},'sss')
        if exist('omega', 'var') && ~isempty(omega)
            % --------- frequency range values given ---------
            if varargin{i}.Ts == 0
                m = freqrespCell(varargin{i},1i* omega);
            else
                m = freqrespCell(varargin{i},exp(1i* omega*varargin{i}.Ts));
            end
        else
            % TODO: fix for Ts~=0
            [m, omega] = getFreqRange(varargin{i});
        end
        
        if isempty(varargin{i}.Name)
            varargin{i}.Name = inputname(i);
        end
        
        % store system in caller workspace
        if inputname(1)
            assignin('caller', inputname(i), varargin{i});
        end
        
        % create frequency response data model
        temp = zeros(varargin{i}.p,varargin{i}.m,length(omega));
        for iO = 1:varargin{i}.p
            for iI = 1:varargin{i}.m
                temp(iO,iI,:) = m{iO,iI};
            end
        end
        
        varargin{i} = frd(temp,omega,varargin{i}.Ts,...
            'InputName',varargin{i}.InputName,'OutputName',varargin{i}.OutputName,...
            'Name',varargin{i}.Name);
        
        % clear omega if input consists of more than one system without omega
        if ~nnz(omegaIndex) && nSys>1
            omega=[];
        end
    end
end

% plot
if nargout== 1
    varargout{1} =  varargin{1};
elseif nargout
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}] = bode(varargin{:});
else
    bode(varargin{:});
end
end

function m = freqrespCell(varargin)
m = freqresp(varargin{1},varargin{2});
warning('off', 'MATLAB:mat2cell:TrailingUnityVectorArgRemoved');
m = mat2cell(m,ones(size(m,1),1),ones(size(m,2),1),size(m,3));
m = cellfun(@(x) x(:,:),m, 'UniformOutput', false);
end

function [m,omega] = getFreqRange(sys)

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
        wmin=wmin-1; t = freqrespCell(sys, 1i*10^wmin);
    end
    while cellfun(@(x,y) norm(x-y)/norm(y),t,dc) < 1e-2
        wmin=wmin+1; t = freqrespCell(sys, 1i*10^wmin);
    end
    wmin=wmin-1;
end

%determine maximum frequency
wmax=0; t = freqrespCell(sys, 10^wmax);
while cellfun(@(x,y,z) norm(x-y)/norm(z),t,ft,dc) > 1e-6
    wmax=wmax+1; t = freqrespCell(sys, 1i*10^wmax);
end
while cellfun(@(x,y,z) norm(x-y)/norm(z),t,ft,dc) < 1e-6
    wmax=wmax-1; t = freqrespCell(sys, 1i*10^wmax);
end
wmax=wmax+1;

delta = (wmax-wmin)/19; % initial resolution (insert odd number only!)
omega = 10.^(wmin:delta:wmax);
m = freqrespCell(sys, 1i* omega);

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
    temp=freqrespCell(sys, 1j*omega(2:2:length(omega)));
    
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
