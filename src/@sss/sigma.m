function  [mag, omega] = sigma(sys, varargin)
% sigma - Plots the amplitude of the frequency response of an LTI system
% 
% Syntax: 
%       [mag, omega] = sigma(sys, omega)
%       [mag, omega] = sigma(sys, omega, options)
%
% Description:
%       Plots the amplitude of the frequency response of an LTI system
%
% Input Arguments:       
%       *Required Input Arguments*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments*
%       -omega:     vector of imaginary frequencies to plot at
%       -options:   plot options. see <a href="matlab:help plot">PLOT</a>
%
% Output Arguments:      
%       -mag:   vector of complex frequency response values
%       -omega: vector of complex frequencies
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
% Authors:      Heiko Panzer, Sylvia Cremer, Rudy Eid
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

options = {};
% in=[]; out=[];

% --------------- EVALUATE OPTIONS ---------------

if nargin>1
    options = varargin;
    if strcmp(class(varargin{1}), 'double')
        omega = varargin{1};
        options(1) = [];
    end
%     if length(options)>=2
%         if strcmp(class(options{1}), 'double') && ...
%            strcmp(class(options{2}), 'double')
%             % I/O
%             if ~isempty(options{1}) && ~isnan(options{1}) 
%                 %input number
%                 in=options{1};
%             end
%             if ~isempty(options{2}) && ~isnan(options{2})
%                 %output number
%                 out=options{2};
%             end
%             options(1:2) = [];
%         end
%     end
end

% --------------- CALCULATE FREQUENCY RESPONSE ---------------


if exist('omega', 'var') && ~isempty(omega)
    % --------- frequency values given ---------
    if size(omega,1)>1
        omega=transpose(omega);
    end
    m = freqrespCell(sys, abs(omega));
else
    % --------- frequency values need to be chosen ---------
    dc = freqrespCell(sys,0);    % G(0)=DCgain
    ft = freqrespCell(sys,inf);  % G(inf)=feedthrough
    mx = cellfun(@(x,y,z) max([norm(x),norm(y),norm(z),1e-8]), dc, ft, ...
        freqrespCell(sys,1), 'UniformOutput',0);

    
    wmax = round(log10(abs(eigs(sys.A,sys.E,1,'LM',struct('tol', 1e-4)))));
    
    %determine maximum frequency
    t = freqrespCell(sys, -1i*10^wmax);
    while cellfun(@(x,y,z) min([norm(x-y) norm(x)])/norm(z),t,ft,mx) > 1e-3
        wmax=wmax+1; t = freqrespCell(sys, 10^wmax);
        mx = cellfun(@(x,y) max([norm(x),norm(y)]), t, mx,'UniformOutput',0);
    end
    while cellfun(@(x,y,z) min([norm(x-y) norm(x)])/norm(z),t,ft,mx) < 1e-3
        wmax=wmax-1; t = freqrespCell(sys, 10^wmax);
        mx = cellfun(@(x,y) max([norm(x),norm(y)]), t, mx,'UniformOutput',0);
    end
    wmax=wmax+1;
    
    
    %determine minimum frequency
    if any(any(cellfun(@isinf,dc))) || any(any(cellfun(@isnan,dc))) % pole at s=0
        wmin = 0;   %***
%         dc = num2cell(ones(size(dc)));
    elseif any(any(cellfun(@abs,dc)<1e-10))   % transfer zero at s=0
        wmin = round(log10(abs(eigs(sys.A,sys.E,1,'SM',struct('tol', 1e-4)))));
        t = freqrespCell(sys, 10^wmin);  %***
        while cellfun(@(x,y,z) min([norm(x-y) norm(x)])/norm(z),t,dc,mx) > 1e-4
            wmin=wmin-1; t = freqrespCell(sys, 10^wmin);
            mx = cellfun(@(x,y) max([norm(x),norm(y)]), t, mx,'UniformOutput',0);
        end
        while cellfun(@(x,y,z) min([norm(x-y) norm(x)])/norm(z),t,dc,mx) < 1e-4
            wmin=wmin+1; t = freqrespCell(sys, 10^wmin);
            mx = cellfun(@(x,y) max([norm(x),norm(y)]), t, mx,'UniformOutput',0);
        end
%         wmin=wmin-1;
%         dc = num2cell(zeros(size(dc)));
    else
        wmin = round(log10(abs(eigs(sys.A,sys.E,1,'SM',struct('tol', 1e-4)))));
        t = freqrespCell(sys, -1i*10^wmin);
        while cellfun(@(x,y,z) min([norm(x-y) norm(x)])/norm(z),t,dc,mx) > 1e-2
            wmin=wmin-1; t = freqrespCell(sys, 10^wmin);
            mx = cellfun(@(x,y) max([norm(x),norm(y)]), t, mx,'UniformOutput',0);
        end
        while cellfun(@(x,y,z) min([norm(x-y) norm(x)])/norm(z),t,dc,mx) < 1e-2
            wmin=wmin+1; t = freqrespCell(sys, 10^wmin);
            mx = cellfun(@(x,y) max([norm(x),norm(y)]), t, mx,'UniformOutput',0);
        end
        wmin=wmin-1;
    end
    

    delta = (wmax-wmin)/39; % initial resolution (insert odd number only!)
    omega = 10.^(wmin:delta:wmax);
    m = freqrespCell(sys, omega);
    
    idx = 2:(length(omega)-1);
    while(1)
        % increase plot density until vertical change per step is below 1%
        newomega = [];
        for i=idx
            if cellfun(@(x) abs( 2*log(abs(x(i))) - log(abs(x(i-1))) - log(abs(x(i+1))) )/(log(max(x))-log(min(x))),m) > 0.001
                newomega = [newomega (omega(i-1)+omega(i))/2 (omega(i)+omega(i+1))/2];
%                 newomega(i-1)=NaN;
            end
        end
        if isempty(newomega), break, end

        % calculate new values of frequency response
        temp=freqrespCell(sys, newomega);
        
%         loglog(abs(omega),abs(m{:}), '.')
%         hold on
%         loglog(abs(newomega),abs(temp{:}), 'r.')
        
        [omega,idx] = sort([omega newomega]);

        % merge cell arrays
        m=cellfun(@(x,y) reshape([x(:);y(:)],1,length(x)+length(y)),m,temp,'UniformOutput',false);
        % ... and sort them
        m=cellfun(@(x) x(idx), m, 'UniformOutput',false);
        
        idx = find(idx > length(omega)-length(newomega));
        
        % do not refine above 1000 points
        if length(omega)>1000
            break
        end
    end
end

% determine magnitude and phase from complex frequency response
mag = cellfun(@abs, m, 'UniformOutput',false);
% phase = cellfun(@(x) angle(x)*180/pi,m, 'UniformOutput',false);

% determine H_inf-norm (maximum magnitude)
[a,b]=cellfun(@max, mag);
[a,c]=max(a);
[H_inf_norm,d]=max(a);
H_inf_peakfreq=omega(b(c(d),d));
if any(H_inf_norm > sys.H_inf_norm) || isempty(sys.H_inf_norm)
    sys.H_inf_norm = H_inf_norm;
    sys.H_inf_peakfreq = H_inf_peakfreq;
end
% store system in caller workspace
if inputname(1)
    assignin('caller', inputname(1), sys);
end

% --------------- PLOT ---------------
if nargout>0
    return
end

if isempty(options)
    options={'Color','b'};
end

mag=cellfun(@(x) 20*log10(x), mag, 'UniformOutput',false);
[m, p] = size(mag); % inputs, outputs
plot_handles=zeros(m,p);

for i=1:m %for each input
    for j=1:p %for each output
        plot_handles(i,j)=subplot(m,p,((i-1)*p)+j);
        hold on
        if j==1 && sys.p>1
            y_lab=sprintf('To Out(%i)',ceil(i/2));
            ylabel(y_lab)
        end
        if i==1 && sys.m>1
            x_lab=sprintf('From In(%i)',j);
            title(x_lab)
        end
      
        % amplitude
        y_plot=mag{i,j};
        plot(omega,y_plot,options{:})
        set(gca, 'XScale', 'log');
        set(gca, 'XLim', [min(omega) max(omega)]);
        mx=max(mag{i,j}); mn=min(mag{i,j});
        if isinf(mx)
            warning('System is zero.');
            continue
        end
        set(gca, 'YLim', [mn-(mx-mn)/20, mx+(mx-mn)/20]);
        ylabel('Magnitude [dB]') 
        xlabel('Frequency [rad/sec]') 
    end
end

% figure and axes properties
set(gcf,'UserData',plot_handles)

h= plot_handles;
% h=axes('position',[0,0,1,1],'Visible','off');
% text(0.4,0.98,'Bode Diagram');
% text(0.5,0.02,'Frequency [rad/sec]')
% text(0.01,0.5,'Magnitude [dB] ; Phase [deg]','Rotation',90)
% set(h,'HandleVisibility','off')

% for i=1:size(plot_handles,1)
%     for j=1:size(plot_handles,2)
%         axes(plot_handles(i,j))
%     end
% end
% axes(plot_handles(1,1))

if nargout==0
    clear mag phase w
end
end

function [m, omega] = freqrespCell(varargin)
    [m, omega] = freqresp(varargin{:});
    warning('off', 'MATLAB:mat2cell:TrailingUnityVectorArgRemoved');
    m = mat2cell(m,ones(size(m,1),1),ones(size(m,2),1),size(m,3));
    m = cellfun(@(x) x(:,:),m, 'UniformOutput', false);
end