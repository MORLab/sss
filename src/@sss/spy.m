function spy(sys,name)
% SPY - Plot sparsity pattern of sss system
% 
% Syntax:
%       SPY(sys)
%       SPY(sys,name)
% 
% Description:
%       This function plots the sparsity pattern of the E and A matrices of
%       the sparse state-space system sys. 
%
%       It is possible to pass a name via a second optional argument.
%
% Input Arguments:
%       *Required Input Arguments:*
%		-sys:  sparse state space (sss)-object
%       *Optional Input Arguments:*
%       -name: Plot title
%       
% Output Arguments:
%       No Output Arguments, just a figure
%
% Examples:
%		The following code plots the sparsity pattern of the benchmark
%		'SpiralInductorPeec' (DSSS, SISO):
%
%> sys = loadSss('SpiralInductorPeec.mat'); 
%> spy(sys,'Peec inductor');
%
%//Note: The .mat file for the example can be found in the benchmarks folder
%
% See Also: 
%		spy
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

    figure;
    subplot(1,2,1);spy(sys.E); title('spy(E)');
    subplot(1,2,2);spy(sys.A); title('spy(A)');
    
    if nargin > 1
        onetitle(name);
    elseif ~isempty(sys.Name)
        onetitle(sys.Name);
    end
end

function onetitle(str)
    %   Create one common title for different subplots
    set(gcf,'NextPlot','add');
    ha = axes; h = title(str);
    set(ha,'Visible','off');
    set(h,'Visible','on'); 
end