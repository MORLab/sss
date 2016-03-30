function  [s, omega] = sigma(sys, varargin)
% sigma - Plots the amplitude of the frequency response of an LTI system
% 
% Syntax: 
%       [s, omega] = sigma(sys, omega)
%       [s, omega] = sigma(sys, omega, options)
%
% Description:
%       Plots the singular values of the frequency response of an LTI system
%
% Input Arguments:       
%       *Required Input Arguments*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments*
%       -omega:     vector of imaginary frequencies to plot at
%       -options:   plot options. see <a href="matlab:help plot">PLOT</a>
%
% Output Arguments:      
%       -s:     vector of singular values of complex frequency response
%       -omega: vector of complex frequencies
%
% Examples:
%       The following code plots the singular values of the frequency 
%       response of the benchmark 'CDplayer' (SSS, MIMO):
%
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> sigma(sys);
%
% See Also:
%       bode, freqresp, magPlot
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
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
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
s=zeros(size(mag,1),size(mag,3));
for i=1:length(omega)
    s(:,i)=svd(mag(:,:,i));
end

if nargout>0 %no plot
    return
end

%% Plot
if isempty(options)
    options={'Color','b'};
end

s=mag2db(s);
for i=1:sys.m
    semilogx(omega,s(i,:),options{:})
    hold on;
end
set(gca, 'XLim', [min(omega) max(omega)]);
xlabel('Frequency (rad/s)');
ylabel('Singular Values (dB)')
title('Singular Values');
end
