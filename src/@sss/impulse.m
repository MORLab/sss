function  [temp,t] = impulse(sys, varargin)
% IMPULSE - Computes and/or plots the impulse response of a sparse LTI system
%
% Syntax:
%       [h, t] = impulse(sys,t)
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
%       The following code computes the impulse response of the benchmark
%       'building' (SSS, SISO) and compares it with the MATLAB built-in function:
%
%> load building.mat
%> sysSparse=sss(A,B,C); %sparse state-space (sss)
%> sys=ss(sysSparse); %full state-space (ss)
%> figure; impulse(sys); hold on; impulse(sysSparse);
%> legend('ss/impulse','sss/impulse');
%
% See Also:
%       residue, ss/impulse, step
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
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
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
if nargin>1 && isa(options{1},'double')
    t=varargin{1};
    options(1)=[];
end


if builtinMATLAB %running the MATLAB built-in function step
    if exist('t','var')
        if nargout>0
            [temp,t]=impulse(ss(sys),t);
        else
            impulse(ss(sys),t);
        end
    else
        if nargout>0
            [temp,t]=impulse(ss(sys),t);
        else
            impulse(ss(sys),t);
        end
    end
else %compute the impulse response with the help of the residue function
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
        % Change the format of res to work with the following code
        
        if size(t,1)>1 && size(t,2)==1
            t=t';
        elseif size(t,1)>1 && size(t,2)>1
            error('t must be a vector.');
        end
        
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
        if isnan(tmax)||isinf(tmax)
            tmax=100;
        end
        delta = tmax/999;
        t = 0:delta:tmax;
        
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
        temp=zeros(length(t),sys.p,sys.m);
        for i=1:sys.p
            for j=1:sys.m
                temp(:,i,j)=h{i,j}';
            end
        end
        
        
        %varargout{1}=temp;
        %varargout{2}=t;
        t=t';
        return
    end
    
    % --------------- PLOT ---------------
    %Defining axis limits
    maxOutput=max(cellfun(@max,h),[],2);
    minOutput=min(cellfun(@min,h),[],2);
    deltaOutput=0.1*(maxOutput-minOutput);
    orderMagnitude=floor(log10(deltaOutput))-1;
    
    heightAxis=round(deltaOutput.*10.^(-orderMagnitude)).*10.^orderMagnitude;
    minOutputAxis=minOutput-heightAxis/2;
    maxOutputAxis=maxOutput+heightAxis/2;
    
    minOutputAxis(minOutput>0&minOutputAxis<=0)=0;
    maxOutputAxis(maxOutput<0&maxOutputAxis>=0)=0;
    %Generating figure to plot
    graphic=stepplot(ss(zeros(1,1),zeros(1,sys.m),zeros(sys.p,1),zeros(sys.p,sys.m)),[-2,-1]);
    opt=getoptions(graphic);
    opt.Title.String='Impulse Response';
    setoptions(graphic,opt);
    
    fig_handle=gcf;
    % set random color if figure is not empty
    if isempty(options)
        if ~isempty(get(fig_handle, 'Children'))
            c=rand(3,1); c=c/norm(c);
            options = {'Color', c};
        end
    end
    
    for jIn=1:sys.m
        for iOut=1:sys.p
            ax=fig_handle.Children(sys.p*sys.m+1-(jIn-1)*sys.p-(iOut-1));
            ax.NextPlot='add';
            
            plot(ax,t, h{iOut,jIn}, options{:});
            set(ax, 'XLim', [0 max(t)], 'YLim', [minOutputAxis(iOut),maxOutputAxis(iOut)]);
        end
    end
    % avoid output
    clear temp h t
    
end
