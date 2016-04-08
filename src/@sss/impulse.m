function  varargout = impulse(varargin)
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
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  10 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% final Time
t = [];
tIndex = cellfun(@isfloat,varargin);
if ~isempty(tIndex) && nnz(tIndex)
    t = varargin{tIndex};
    varargin(tIndex)=[];
end

Tfinal = 0;
for i = 1:length(varargin)
    % Set name to input variable name if not specified
    if isprop(varargin{i},'Name')
        if isempty(varargin{i}.Name) % Cascaded if is necessary && does not work
            varargin{i}.Name = inputname(i);
        end
    end
    
    % Convert sss to frequency response data model
    if isa(varargin{i},'sss')
        [varargin{i},t] = gettf(varargin{i}, t);
    end
    Tfinal = max(t(end),Tfinal);
end

% Call ss/impulse
if nargout == 1 && Opts.frd
    varargout{1} = varargin{1};
elseif nargout
    [varargout{1},varargout{2},varargout{3},varargout{4}] = impulse(varargin{:},Tfinal);
else
    impulse(varargin{:},Tfinal);
end

function [TF,ti] = gettf(sys, t)

[h,t] = impulseLocal(sys, t);
ti = linspace(t(1),t(end),length(t));

Ts = min(diff(ti));
h_ = cell([size(h,2) size(h,3)]);
for i = 1:size(h,2)
    for j = 1:size(h,3)
        h_{i,j} = interp1(t,h(:,i,j),ti)*Ts;        
    end
end

TF = filt(h_,1,Ts,...
    'InputName',sys.InputName,'OutputName',sys.OutputName,...
    'Name',sys.Name);

function [temp,t] = impulseLocal(sys, t)

[res,p]=residue(sys);
builtinMATLAB=1;
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

if builtinMATLAB %running the MATLAB built-in function step            
    [temp,t]=impulse(ss(sys),t);
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
    if ~isempty(t) && length(t)>1
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
        
        if isempty(t)
            tmax = decayTime(sys);
        else
            tmax = t;
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
            if ex==2
                break
            end
        end
    end
    
    temp=zeros(length(t),sys.p,sys.m);
    for i=1:sys.p
        for j=1:sys.m
            temp(:,i,j)=h{i,j}';
        end
    end
    t=t';
    
end
