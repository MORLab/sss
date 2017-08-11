function fNew = bodeToFig(closeBode,defPos,subM,subN)
% BODETOFIG - Transforms MATLAB bode/step/impulse plot into normal figure
%
% Syntax: 
%       BODETOFIG(closeBode,defPos,subM,subN)
%
% Description:
%       Transforms MATLAB bode/step/impulse plot into normal figure
%
% Input Arguments:
%           -closeBode: closes bode figure if true
%           -defPos:    sets axis position to default if true
%           -subM:      m as in subplot
%           -subN:      n as in subplot
%
% Output Arguments:
%           -fNew:     handle to figure
%
% See Also:
%       
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Professur fuer Thermofluiddynamik, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Stefan Jaensch (Jaensch@tfd.mw.tum.de)
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  11 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if ~nargin
    closeBode = false; defPos = false;
end

fig = gcf;
f_new = figure;
for i = 1:length(fig.Children)
    childFig = fig.Children(i);
    switch childFig.Type
        case 'axes'
            %copy only visible axes
            if strcmp(childFig.Visible,'on')
                ax_new = copyobj(childFig,f_new);
                if defPos
                    set(ax_new,'Position','default')
                end
            end
    end
end

if nargin==4
    warning('compare bode and resulting figure')        
    k = 1;
    for i = 1:length(f_new.Children)
        y = f_new.Children(i).Position(2);
        for j = 1:length(f_new.Children)
            if k > length(f_new.Children)
                break
            elseif y == f_new.Children(j).Position(2)
                subplot(subM,subN,k,f_new.Children(j))
                k = k+1;
            end
        end
    end
end

if closeBode
    close(fig);
end