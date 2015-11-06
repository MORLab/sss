function [p, z] = pzmap(sys, varargin)
% PZMAP - Pole-zero plot of sparse state-space system
%
% Syntax:
%       PZMAP(sys)
%       PZMAP(sys, opts)
%       [p,z] = PZMAP(sys)
%
% Description:
%       pzmap(sys) creates a pole-zero plot of the continuous- or discrete-time 
%       dynamic system model sys. For MIMO systems, pzmap plots the system poles
%       and invariant zeros of each SISO model in an individual subplot. The poles 
%       are plotted as x's and the invariant zeros are plotted as o's.
%
%       [p,z] = pzmap(sys) returns the system poles and invariant zeros in the column 
%       vectors p and z. No plot is drawn on the screen.
%
% Input Arguments:
%       -sys:      an sss-object containing the LTI system
%       -opts:     plot options. see <a href="matlab:help plot">PLOT</a>
%
% Output Arguments:
%       -p: vector containing poles 
%       -z: vector containing invariant zeros
%
% Examples:
%       Create a random descriptor model (DSS, SISO) and compare the output
%       of ss/pzmap and sss/pzmap:
%
%> A = randn(500,500); B = randn(500,1); C = randn(1,500); D = zeros(1,1);
%> E = randn(500,500);
%> sys = dss(A,B,C,D,E);
%> sysSss = sss(sys);
%> figure; pzmap(sys);
%> figure; pzmap(sysSss);
%
%       Load the benchmark "PEEC_MTLn1600" (DSS,MIMO) and use pzmap:
%
%> load PEEC_MTLn1600.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m),E)
%> sysTrunc = sys(1:3,1:2);
%> figure; pzmap(sysTrunc);
%
% See Also:
%       ss/pzmap
%
% References:
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Sylvia Cremer, Maria Cruz Varona
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

% are poles already available?
if ~isempty(sys.poles)
    p=sys.poles;
else
    % no, they are not. solve for eigenvalues of system
    p = eig(sys);
    
    sys.poles=p;
    
    % store system in caller workspace
    if inputname(1)
        assignin('caller', inputname(1), sys);
    end
end
% ensure column vector
if size(p,1)<size(p,2)
    p=transpose(p);
end

% are zeros already available?
if ~isempty(sys.invariantZeros)
    z=sys.invariantZeros;
else
    % no, they are not. solve for gen. eigenvalues of Rosenbrock matrix
    %z = cell(sys.p,sys.m);
    %for i=1:sys.m
    %    for j=1:sys.p
    if sys.m==sys.p
            z=eig(full([sys.A,sys.B;sys.C,sys.D]), ...
                       [full(sys.E),zeros(sys.n,sys.m);zeros(sys.p,sys.n),zeros(sys.p,sys.m)]);
    else
        z=zeros(0,1);
    end
            % ensure column vector
     %       if size(z{j,i},1)<size(z{j,i},2)
     %           z{j,i}z=transpose(z{j,i});
     %       end
     %   end
    %end

    % remove zeros at infinity
    %z=cellfun(@(x) x(~isinf(x)), z, 'UniformOutput', false);
    z=z(~isinf(z));
    sys.invariantZeros=z;
    
    % store system in caller workspace
    if inputname(1)
        assignin('caller', inputname(1), sys);
    end
end

if nargout>0
    return
end

% --------------- PLOT ---------------

options=varargin;
fig_handle=gcf; %generate figure
%axes_handle=zeros(sys.p,sys.m);

% set random color if figure is not empty
if isempty(options)
    if ~isempty(get(fig_handle, 'Children'))
        c=rand(3,1); c=c/norm(c);
        options = {'Color', c};
    end
end
% loop for plotting the pole-zero-maps
%for j=1:sys.p %secondly, go through all outputs
%    for i=1:sys.m %firstly, go through all inputs
        %rows: outputs, columns: inputs    
        %axes_handle(j,i)=subplot(sys.p,sys.m,j*sys.m+i-sys.m);
        axes_handle=subplot(1,1,1);
        box on
        
        % plot o for zeros
        %plot_handle=plot(real(z{j,i}), imag(z{j,i}), options{:});
        plot_handle=plot(real(z), imag(z), options{:}); 
        set(plot_handle, 'Marker', 'o', 'LineStyle', 'none');
        hold on
        % plot x for poles, remove legend entry
        plot_handle=plot(real(p), imag(p), options{:});
        set(plot_handle, 'Marker', 'x', 'LineStyle', 'none');
        hAnnotation = get(plot_handle,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')

        % determine x and y boundary
      %  mni=min(imag([p;z{j,i}])); mxi=max(imag([p;z{j,i}]));
        mni=min(imag([p;z])); mxi=max(imag([p;z]));
        if mni*mxi<0 
            limy = [mni-(mxi-mni)/20 mxi+(mxi-mni)/20];
        elseif mni*mxi>0 
            limy = sort([0 1.05*max(abs([mni mxi]))]*sign(mni));
        end
        if sys.Ts ~= 0 %adjust the y-axis so that the unitary circle is visible
            if limy(1)>-1
                limy(1) = -1;
            end
            if limy(2)<1
                limy(2) = 1;
            end
        end
        
        %mnr=min(real([p;z{j,i}])); mxr=max(real([p;z{j,i}]));
        mnr=min(real([p;z])); mxr=max(real([p;z]));
        limx=[mnr*1.05,mxr*1.05];
        if mnr*mxr<0 
            limx = [mnr-(mxr-mnr)/20 mxr+(mxr-mnr)/20];
        elseif mnr*mxr>0
            limx = sort([0 1.05*max(abs([mnr mxr]))]*sign(mnr));
        end
        if sys.Ts ~= 0 %adjust the x-axis so that the unitary circle is visible
            if limx(1)>-1
                limx(1) = -1;
            end
            if limx(2)<1
                limx(2) = 1;
            end
        end
        % plot dashed axes through origin, remove legend entry
        plot_handle=plot([0 0],limy,':k');
        hAnnotation = get(plot_handle,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')
        plot_handle=plot(limx,[0 0],':k');
        hAnnotation = get(plot_handle,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')
        set(axes_handle(1,1), 'XLim', limx, 'YLim', limy);
        
        if sys.Ts ~= 0 %plot unitary circle in case of discrete system
            r=1; %radius
            circlePoints=0:0.01:2*pi;
            plot_handle=plot(r*cos(circlePoints),r*sin(circlePoints),':k');
        end

        % label input / output number
%         if j==1 && sys.m>1
%             x_lab=sprintf('From In(%i)',i);
%             title(x_lab)
%         end
%         if i==1 && sys.p>1
%             y_lab=sprintf('To Out(%i)',j);
%             ylabel(y_lab)
%         end
%    end
%end

% create invisible background plot for labelling
% h=axes('position',[0,0,1,1],'Visible','off');
% text(0.4,0.98,'Pole-Zero Map');
% text(0.5,0.02,'Real Axis')
% text(0.01,0.5,'Imaginary Axis','Rotation',90)
% set(h,'HandleVisibility','off')
% hold on
% % bring subplots to front
% for i=1:size(axes_handle,1)
%     for j=1:size(axes_handle,2)
% %         set(fig_handle, 'CurrentAxes', axes_handle(i,j));
%         axes(axes_handle(i,j))
%     end
% end
% set(fig_handle, 'CurrentAxes', axes_handle(1,1));

% make subplots content of the figure
set(fig_handle,'UserData',axes_handle)

% avoid output
clear p z
