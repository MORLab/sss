function  [temp, t] = step(sys, varargin)
% STEP - Computes and/or plots the step response of a sparse LTI system
%
% Syntax:
%       STEP(sys)
%       [h, t] = STEP(sys)
%       [h, t] = STEP(sys, t)
%       [h, t] = STEP(sys, t, opts)
%
% Description:
%       step(sys) plots the step response of the sparse LTI system sys
%
%       [h, t] = step(sys, t) computes the step response of the sparse LTI
%       system sys and returns the vectors h and t with the response and
%       the time series, respectively.
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
%       The following code plots the step response of the benchmark 
%       'CDplayer' (SSS, MIMO):
%
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> step(sys);
%
% See Also:
%       residue, impulse
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Professur fuer Thermofluiddynamik, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Sylvia Cremer, Jorge Luiz Moreira Silva
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  10 Nov 2015
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
builtinMATLAB=0;
if condest(sys.E)>1/(100*eps) %Verify if E is singular
    warning('Matrix E is close to singular or badly scaled: running the MATLAB built-in function step');
    builtinMATLAB=1;
end
if not(builtinMATLAB)
    for i=1:numel(p) %Verify if the poles are repeated
        if p(i)~=0
            numberOfRepeatedPoles=sum(abs((p-p(i))/p(i))<100*eps);
        else
            numberOfRepeatedPoles=sum(abs(p-p(i))<100*eps);
        end
        if numberOfRepeatedPoles>1
            builtinMATLAB=1;
            warning('System is not diagnoalizable because of repetition of eigenvalues: running the MATLAB built-in function step');
            break;
        end
    end
end
options=varargin;
if nargin>1 && isa(options{1}, 'double')
    t=varargin{1};
    options(1)=[];
end

if builtinMATLAB
    if exist('t','var')
        if nargout>0
            [temp,t]=step(ss(sys),t);
        else
            step(ss(sys),t);
        end
    else
        if nargout>0
            [temp,t]=step(ss(sys));
        else
            step(ss(sys));
        end
    end
        
else
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
    % is time vector given?
    if exist('t', 'var') && ~isempty(t)
        % calculate step response
        diagCell = cellfun(@(x) (x./p),res,'UniformOutput',false);
        diagCell=cellfun(@nanzero,diagCell,'UniformOutput',false);
        resZero=cellfun(@(x)x(p==0),res,'UniformOutput',false);
        if numel(resZero{1})==0
            resZero=num2cell(zeros(sys.p,sys.m));
        end
        h=cellfun(@(x,y,z) sum(x*transpose(exp(t'*p)-1),1)+y+z*t,diagCell,num2cell(sys.D),resZero,'UniformOutput',false);
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
        if isnan(tmax)||isinf(tmax)
            tmax=100;
        end
        delta = tmax/999;
        t = 0:delta:tmax;
        
        % calculate step response
        diagCell = cellfun(@(x) (x./p),res,'UniformOutput',false);
        diagCell=cellfun(@nanzero,diagCell,'UniformOutput',false);
        resZero=cellfun(@(x)x(p==0),res,'UniformOutput',false);
        if numel(resZero{1})==0
            resZero=num2cell(zeros(sys.p,sys.m));
        end
        h=cellfun(@(x,y,z) sum(x*transpose(exp(t'*p)-1),1)+y+z*t,diagCell,num2cell(sys.D),resZero,'UniformOutput',false);
        h=cellfun(@real,h,'UniformOutput',false);
        
        % increase resolution ias long as rel. step size is too large
        ex=1;
        while 1
            refine=0;
            for iOut=1:sys.p
                for jIn=1:sys.m
                    m=h{iOut,jIn};
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
        temp=zeros(length(t),sys.p,sys.m);
        for i=1:sys.p
            for j=1:sys.m
                temp(:,i,j)=h{i,j}';
            end
        end
        t=t';
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
    deltaOutput=0.1*(maxOutput-minOutput);
    orderMagnitude=floor(log10(deltaOutput))-1;
    
    heightAxis=round(deltaOutput.*10.^(-orderMagnitude)).*10.^orderMagnitude;
    minOutputAxis=minOutput-heightAxis/2;
    minOutputAxis(isnan(minOutputAxis))=-inf;
    maxOutputAxis=maxOutput+heightAxis/2;
    maxOutputAxis(isnan(maxOutputAxis))=inf;
    
    minOutputAxis(minOutput>0&minOutputAxis<=0)=0;
    maxOutputAxis(maxOutput<0&maxOutputAxis>=0)=0;
    
    
    for iOut=1:sys.p
        for jIn=1:sys.m
            axes_handle(iOut,jIn)=subplot(sys.p,sys.m,(iOut-1)*sys.m+jIn);
            hold on;
            
            if jIn==1 && sys.p>1
                y_lab=sprintf('To Out(%i)',ceil(iOut/2));
                ylabel(y_lab,'FontSize',10,'FontName','Helvetica','Color',[0.31,0.31,0.31],...
                    'FontWeight','normal','FontAngle','normal');
            end
            if iOut==1 && sys.m>1
                x_lab=sprintf('From In(%i)',jIn);
                title(x_lab,'FontSize',10,'FontName','Helvetica','Color',[0.31,0.31,0.31],...
                    'FontWeight','normal','FontAngle','normal');
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
    yLabel=text(0.01,0.5,'Amplitude','FontName','Helvetica','FontSize',11,'FontWeight','normal','rotation',90);
    set(yLabel,'HorizontalAlignment','center','VerticalAlignment','top');
    uistack(labelsAxes,'bottom');
    z=zoom;
    setAllowAxesZoom(z,labelsAxes,false);
    clear h t
end
end
function x = nanzero(x)
x(isnan(x)) = 0;
return
end