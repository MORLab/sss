function  [h, t] = step(sys, varargin)
% STEP - Computes and/or plots the step response of a sparse LTI system
% 
% Syntax:
%       STEP(sys)
%       [h, t] = STEP(sys)
%       [h, t] = STEP(sys, t)
%       [h, t] = STEP(sys, t, opts)
%
% Description:
%       Computes and/or plots the step response of a sparse LTI system
% 
% Input Arguments:       
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:* 
%       -t:     vector of time values to plot at
%       -opts:  plot options, see <a href="matlab:help plot">PLOT</a>
%
% Output Arguments:      
%       -h, t: vectors containing step response and time vector
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
% Authors:      Heiko Panzer, Sylvia Cremer, Jorge Luiz Moreira Silva
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
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
    % calculate step response
    h=cellfun(@(x,y) sum(diag(x./p)*transpose(exp(t'*p)-1),1)+y,res,num2cell(sys.D),'UniformOutput',false);
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
    if isnan(tmax)
        tmax=100;
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

    % calculate step response
    diagCell = cellfun(@(x) diag(x./p),res,'UniformOutput',false);
    h=cellfun(@(x,y) sum(x*transpose(exp(t'*p)-1),1)+y,diagCell,num2cell(sys.D),'UniformOutput',false);
    h=cellfun(@real,h,'UniformOutput',false);

    % increase resolution ias long as rel. step size is too large
    ex=1;
    while 1
        refine=0;
        for iOut=1:sys.p
            for jIn=1:sys.m
                m=h{iOut,jIn};
%                 for k=2:length(m)-1
%                     if abs(abs(m(k)) - abs(m(k+1)))/(abs(m(k)) + abs(m(k+1))) > 0.5
                    if any(abs(diff(abs(m)))./((abs(m(1:end-1)))) > 0.5)
                        delta=delta/2;
                        t=0:delta:tmax;
                        t_temp=t(2:2:end);
                        temp=cellfun(@(x,y) sum(x*transpose(exp(t_temp'*p)-1),1)+y,diagCell,num2cell(sys.D),'UniformOutput',false);
                        temp=cellfun(@real,temp,'UniformOutput',false);
                        h=cellfun(@(x,y) [reshape([x(1:length(x)-1); y],1,2*length(x)-2),x(end)],h,temp,'UniformOutput',false);
                        refine=1;
                        break
                    end
%                 end
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

maxOutput=max(cellfun(@max,h),[],2);
minOutput=min(cellfun(@min,h),[],2);
deltaOutput=0.2*(maxOutput-minOutput);
orderMagnitude=floor(log10(deltaOutput))-1;

heightAxis=round(deltaOutput.*10.^(-orderMagnitude)).*10.^orderMagnitude;
minOutputAxis=minOutput-heightAxis/2;
maxOutputAxis=maxOutput+heightAxis/2;

minOutputAxis(minOutput>0&minOutputAxis<=0)=0;
maxOutputAxis(maxOutput<0&maxOutputAxis>=0)=0;


for iOut=1:sys.p
    for jIn=1:sys.m
        axes_handle(iOut,jIn)=subplot(sys.p,sys.m,(iOut-1)*sys.m+jIn);
        hold on;

        if jIn==1 && sys.p>1
            y_lab=sprintf('To Out(%i)',ceil(iOut/2));
            ylabel(y_lab,'FontSize',10,'FontName','Helvetica','Color',[0.31,0.31,0.31],...
            'FontWeight','normal','FontAngle','normal'); %'FontSmoothing','on'
        end
        if iOut==1 && sys.m>1
            x_lab=sprintf('From In(%i)',jIn);
            title(x_lab,'FontSize',10,'FontName','Helvetica','Color',[0.31,0.31,0.31],...
            'FontWeight','normal','FontAngle','normal'); %'FontSmoothing','on',
        end
        if jIn==1 &&(iOut==sys.p)
            %do nothing
        elseif iOut==sys.p
            set(gca,'ytick',[])
        elseif jIn==1
            set(gca,'xtick',[])
        else
            set(gca,'xtick',[],'ytick',[])
        end

        plot(t, h{iOut,jIn}, options{:});
        hold on
        plot([0,max(t)],[h{iOut,jIn}(end),h{iOut,jIn}(end)],':','Color',[0.31 0.31 0.31]);
        axis([0,max(t),minOutputAxis(iOut),maxOutputAxis(iOut)]);
    end
end
labelsAxes=axes('Position',[0 0 1 1],'Visible','off');
text(0.5,.99,'Step Response','FontName','Helvetica','FontSize',11,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','cap');
text(0.5,0.01,'Time (seconds)','FontName','Helvetica','FontSize',11,'FontWeight','normal','HorizontalAlignment','center','VerticalAlignment','bottom');
h=text(0.01,0.5,'Amplitude','FontName','Helvetica','FontSize',11,'FontWeight','normal','rotation',90);
set(h,'HorizontalAlignment','center','VerticalAlignment','top');
uistack(labelsAxes,'bottom');
z=zoom;
setAllowAxesZoom(z,labelsAxes,false);
clear h t