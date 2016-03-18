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
%       dynamic system model sys. For MIMO systems, pzmap plots the system 
%       poles and the invariant zeros in one figure. The poles are plotted 
%       as x's and the invariant zeros are plotted as o's.
%
%       [p,z] = pzmap(sys) returns the system poles and invariant zeros in the column 
%       vectors p and z. No plot is drawn on the screen.
%
%//Note: The calculation of the invariant zeros is only defined for systems
%       with the same number of inputs and outputs (m=p). That means that if
%       pzmap is called with a system with m~=p, then z = [ ].
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
%       Create a random descriptor model (DSSS, SISO) and compare the output
%       of ss/pzmap and sss/pzmap:
%
%> A = randn(500,500); B = randn(500,1); C = randn(1,500); D = zeros(1,1);
%> E = randn(500,500);
%> sys = dss(A,B,C,D,E);
%> sysSss = sss(sys);
%> figure; pzmap(sys);
%> figure; pzmap(sysSss);
%
%       Load the benchmark 'rail_1357' (DSSS, MIMO) and use pzmap:
%
%> load rail_1357.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m),E)
%> figure; pzmap(sys);
%
% See Also:
%       ss/pzmap
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
% Authors:      Heiko Panzer, Sylvia Cremer, Maria Cruz Varona
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  10 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

% no, they are not. solve for eigenvalues of system
p = eig(sys);
    
% ensure column vector
if size(p,1)<size(p,2)
    p=transpose(p);
end


if sys.m==sys.p
        z=eig(full([sys.A,sys.B;sys.C,sys.D]), ...
                   [full(sys.E),zeros(sys.n,sys.m);zeros(sys.p,sys.n),zeros(sys.p,sys.m)]);
else
    z=zeros(0,1);
end

% remove zeros at infinity
z=z(real(-z)<1e11);
z=z(~isinf(z));

if nargout>0
    return
end

% --------------- PLOT ---------------

options=varargin;
fig_handle=gcf; %generate figure

% set random color if figure is not empty
if isempty(options)
    if ~isempty(get(fig_handle, 'Children'))
        c=rand(3,1); c=c/norm(c);
        options = {'Color', c};
    end
end
        axes_handle=subplot(1,1,1);
        box on
        
        % plot o for zeros
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
        mni=min(imag([p;z])); mxi=max(imag([p;z]));
        limy=[mni*1.05,mxi*1.05];
        if mni*mxi<0 
            limy = [mni-(mxi-mni)/20 mxi+(mxi-mni)/20];
        elseif mni*mxi>0 
            limy = sort([0 1.05*max(abs([mni mxi]))]*sign(mni));
        elseif mni==0 && mxi==0
            limy = [-1 1];
        end
        if sys.Ts ~= 0 %adjust the y-axis so that the unitary circle is visible
            if limy(1)>-1
                limy(1) = -1;
            end
            if limy(2)<1
                limy(2) = 1;
            end
        end

        mnr=min(real([p;z])); mxr=max(real([p;z]));
        limx=[mnr*1.05,mxr*1.05];
        if mnr*mxr<0 
            limx = [mnr-(mxr-mnr)/20 mxr+(mxr-mnr)/20];
        elseif mnr*mxr>0
            limx = sort([0 1.05*max(abs([mnr mxr]))]*sign(mnr));
        elseif mnr==0 && mxr==0
            limx = [-1 1];
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
        xlabel('Real Axis (seconds^{-1})');
        ylabel('Imaginary Axis (seconds^{-1}');
        title('Pole-Zero Map');
        
        if sys.Ts ~= 0 %plot unitary circle in case of discrete system
            r=1; %radius
            circlePoints=0:0.01:2*pi;
            plot_handle=plot(r*cos(circlePoints),r*sin(circlePoints),':k');
        end

% make subplots content of the figure
set(fig_handle,'UserData',axes_handle)

% avoid output
clear p z
