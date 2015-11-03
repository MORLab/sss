function [f_new] = bodeToFig(closeBode,defPos,subM,subN)
%bodeToFig - Transforms matlab bode/step/impulse plot into normal figure
%
% Syntax:  bodeToFig(closeBode,defPos,subM,subN)
%
% Inputs:
%    closeBode - closes bode figure if true
%    dePos - sets axis position to default if true
%    subM - m as in subplot
%    subN - n as in subplot,closeBode,defPos
% Outputs:
%   f_new - handle to figure
%
% Other m-files required:
% Subfunctions: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% ------------------------------------------------------------------
% Authors:      Stefan Jaensch (Jaensch@tfd.mw.tum.de)
% ------------------------------------------------------------------

if ~nargin
    closeBode = false; defPos = false;
end

fig = gcf;
f_new = figure;
for i = 1:length(fig.Children)
    childFig = fig.Children(i);
    switch childFig.Type
        case 'axes'
            ax_new = copyobj(childFig,f_new);
            if defPos
                set(ax_new,'Position','default')
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
                subplot( subM,subN,k,f_new.Children(j))
                k = k+1;
            end
        end
    end
end

if closeBode
    close(fig);
end