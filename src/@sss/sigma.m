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
%       The following code plots the amplitude of the frequency response of
%       the benchmark 'CDplayer' (SSS, MIMO):
%
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> sigma(sys);
%
% See Also:
%       bode, freqresp
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
% Authors:      Heiko Panzer, Sylvia Cremer, Rudy Eid, 
%               Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  07 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

options = {};

%% Evaluate options
if nargin>1, options = varargin;
    if isa(varargin{1}, 'double')
        omega = varargin{1}; options(1) = [];
    end
end
%% Calculate frequency response
if exist('omega', 'var') && ~isempty(omega)
    % frequency values given
    m = freqresp(sys, omega);
else
    % frequency values need to be chosen
    [m, omega] = freqresp(sys);
end
mag = abs(m);

if nargout>0 %no plot
    return
end

%% Plot
if isempty(options)
    options={'Color','b'};
end

mag=mag2db(mag); %can be overwritten since it is not returned
m = sys.m; p = sys.p;
plot_handles=zeros(m,p);

for i=1:p %for each output
    for j=1:m %for each input
        plot_handles(i,j)=subplot(p,m,(i-1)+j);
        hold on
        if j==1 && sys.m>1
            y_lab=sprintf('To Out(%i)',ceil(i/2));
            ylabel(y_lab)
        end
        if i==1 && sys.p>1
            x_lab=sprintf('From In(%i)',j);
            title(x_lab)
        end
      
        % amplitude
        y_plot=squeeze(mag(i,j,:));
        plot(omega,y_plot,options{:})
        set(gca, 'XScale', 'log');
        set(gca, 'XLim', [min(omega) max(omega)]);
        mx=max(y_plot); mn=min(y_plot);
        if isinf(mx)
            warning('System is zero.');
            continue
        end
        set(gca, 'YLim', [mn-(mx-mn)/20, mx+(mx-mn)/20]);
        ylabel('Magnitude [dB]') 
        xlabel('Frequency [rad/sec]') 
    end
end

end
