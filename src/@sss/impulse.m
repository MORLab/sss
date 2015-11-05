function  [h, t] = impulse(sys, varargin)
% impulse - Computes and/or plots the impulse response of a sparse LTI system
%
% Syntax:
%       [h, t] = impulse(sys,t,opts)
%
% Description:
%       Computes and/or plots the impulse response of a sparse LTI system
%
% Input Arguments:      
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:* 
%       -t:     vector of time values to plot at
%       -opts:  plot options. see <a href="matlab:help plot">PLOT</a>
%
% Outputs:      
%       -h, t: vectors containing impulse response and time vector
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
% Authors:      Heiko Panzer, Sylvia Cremer
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% are poles and residues already available?
if ~isempty(sys.poles) && ~isempty(sys.residues)
    p=sys.poles;
    res=sys.residues;
else
    [res,p]=residue(sys);
    % store system to caller workspace
    if inputname(1)
        assignin('caller', inputname(1), sys);
    end
end

options=varargin;
if nargin>1 && strcmp(class(options{1}), 'double')
    t=varargin{1};
    options(1)=[];
end

% is time vector given?
if exist('t', 'var') && ~isempty(t)
    % calculate impulse response
    h=cellfun(@(x) sum((diag(x)*conj(exp(t'*p))'),1),res,'UniformOutput',false);
    h=cellfun(@real,h,'UniformOutput',false);
    
else
    % no, retrieve time values automatically
    
    % is decayTime already available?
    if ~isempty(sys.decayTime)
        tmax = sys.decayTime;
    else
        tmax = decayTime(sys);
        % store system to caller workspace
        if inputname(1)
            assignin('caller', inputname(1), sys);
        end        
    end
    delta = tmax/999;
    t = 0:delta:tmax;
    
    % Change the format of res to work with the following code
    %   original: res = {res1, res2,... resN} where resK = res(p(k)) is the
    %           (possibly matrix-valued) residual to p(k)
    %   new     : res = cell(sys.p,sys.m), where res{i,j} is a vector of
    %           scalar residual: res{i,j} = [res{i,j}1, ..., res{i,j}N]
    resOld = res; clear res; res = cell(sys.p,sys.m);
    for iOut = 1:sys.p
        for jIn = 1:sys.m
            res{iOut,jIn} = [];
            for kPole = 1:length(p)
                res{iOut,jIn} = [res{iOut,jIn}, resOld{kPole}(iOut,jIn)];
            end
        end
    end

    % calculate impulse response
%     cellfun(@(x) sum((diag(x)*conj(exp(t'*p))'),1),res,'UniformOutput',false);
    h=cellfun(@(x)   sum(diag(x)*transpose(exp(t'*p)),1),res,'UniformOutput',false);
    h=cellfun(@real,h,'UniformOutput',false);

    % increase resolution as long as rel. step size is too large
    ex=1;
    while 1
        refine=0;
        for iOut=1:sys.p
            for jIn=1:sys.m
                m=h{iOut,jIn};
                for k=2:length(m)-1
                    if abs(abs(m(k)) - abs(m(k+1)))/(abs(m(k)) + abs(m(k+1))) > 0.5
                        delta=delta/2;
                        t=0:delta:tmax;
                        t_temp=t(2:2:end);
                        temp=cellfun(@(x) sum((diag(x)*conj(exp(t_temp'*p))'),1),res,'UniformOutput',false);
                        temp=cellfun(@real,temp,'UniformOutput',false);
                        h=cellfun(@(x,y) [reshape([x(1:length(x)-1); y],1,2*length(x)-2),x(end)],h,temp,'UniformOutput',false);
                        refine=1;
                        break
                    end
                end
                if refine
                    break
                end
            end
            if refine
                break
            end
        end
        if ~refine
            break
        end
        ex=ex+1;
        if ex==5 
            break
        end
    end
end

if nargout>0
    return
end

% --------------- PLOT ---------------

% set random color if figure is not empty
fig_handle=gcf;
if isempty(options)
    if ~isempty(get(fig_handle, 'Children'))
        c=rand(3,1); c=c/norm(c);
        options = {'Color', c};
    end
end

axes_handle=zeros(sys.p,sys.m);
for iOut=1:sys.p
    for jIn=1:sys.m
        axes_handle(iOut,jIn)=subplot(sys.p,sys.m,(jIn-1)*sys.p+iOut);
        hold on;
        
        if jIn==1 && sys.p>1
            y_lab=sprintf('To Out(%i)',ceil(iOut/2));
            ylabel(y_lab,'FontSize',10,'FontName','Helvetica','Color',[0.31,0.31,0.31],...
            'FontWeight','normal','FontSmoothing','on','FontAngle','normal');
        end
        if iOut==1 && sys.m>1
            x_lab=sprintf('From In(%i)',jIn);
            title(x_lab,'FontSize',10,'FontName','Helvetica','Color',[0.31,0.31,0.31],...
            'FontWeight','normal','FontSmoothing','on','FontAngle','normal');
        end
        
        plot(t, h{iOut,jIn}, options{:});
        mx=max(h{iOut,jIn}); mn=min(h{iOut,jIn});
        if mx==0 && mn==0
            mx=1; mn=-1;
        end
        set(gca, 'XLim', [0 max(t)], 'YLim', [mn-(mx-mn)/20, mx+(mx-mn)/20]);

%         ylabel('Amplitude') 
%         xlabel('Time (seconds)') 
        
    end
end

axes('Position',[0 0 1 1],'Visible','off');
text(0.5,.99,'Impulse Response','FontName','Helvetica','FontSize',11,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','cap');
text(0.5,0.01,'Time (seconds)','FontName','Helvetica','FontSize',11,'FontWeight','normal','HorizontalAlignment','center','VerticalAlignment','bottom');
h=text(0.01,0.5,'Amplitude','FontName','Helvetica','FontSize',11,'FontWeight','normal','rotation',90);
set(h,'HorizontalAlignment','center','VerticalAlignment','top');

% avoid output
clear h t


    