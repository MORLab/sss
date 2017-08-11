function sys = append(varargin)
% APPEND - Appends a set of sparse LTI systems (sss)
%
% Syntax:
%       sys = APPEND(sys1,sys2,...)
%
% Description:
%       Create an aggregated model (sss) by concatenating inputs and
%       outputs
%
% Input Arguments:
%       -sys1,sys2,...: sparse state space (sss)-objects
%
% Output Arguments:
%       -sys: appended (open loop) sparse state-space (sss)-object
%
% Examples:
%		This code loads two benchmark models included in the toolbox
%		and groups the models by appending their inputs and outputs:
%
%> load building; 
%> sysBuilding = sss(A,B,C)
%> load CDplayer
%> sysCdplayer = sss(A,B,C)
%> sysAppended = append(sysBuilding,sysCdplayer)
%
% See Also:
%       connect, ss/append, ss/blkdiag
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
% Authors:      Thomas Emmert (emmert@tfd.mw.tum.de)
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  04 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

sys = varargin{1};

for i = 2: nargin
    % Merge input and output groups
    InGroups = varargin{i}.InputGroup;
    if not(isempty(InGroups))
        for Group = fieldnames(InGroups)'
            if not(isempty(Group))
                if isfield(sys.InputGroup,Group)
                    sys.InputGroup.(char(Group)) =  [sys.InputGroup.(char(Group)) ,sys.m + InGroups.(char(Group))];
                else
                    sys.InputGroup.(char(Group)) =  sys.m + InGroups.(char(Group));
                end
            end
        end
    end
    OutGroups = varargin{i}.OutputGroup;
    if not(isempty(OutGroups))
        for Group = fieldnames(OutGroups)'
            if not(isempty(Group))
                if isfield(sys.OutputGroup,Group)
                    sys.OutputGroup.(char(Group)) =  [sys.OutputGroup.(char(Group)) ,sys.p + OutGroups.(char(Group))];
                else
                    sys.OutputGroup.(char(Group)) =  sys.p + OutGroups.(char(Group));
                end
            end
        end
    end
    StateGroups = varargin{i}.StateGroup;
    if not(isempty(StateGroups))
        for Group = fieldnames(StateGroups)'
            if not(isempty(Group))
                if isfield(sys.StateGroup,Group)
                    sys.StateGroup.(char(Group)) =  [sys.StateGroup.(char(Group)) ,sys.n + StateGroups.(char(Group))];
                else
                    sys.StateGroup.(char(Group)) =  sys.n + StateGroups.(char(Group));
                end
            end
        end
    end

    y = [sys.y;varargin{i}.y];
    u = [sys.u;varargin{i}.u];
    StateName = [sys.StateName; varargin{i}.StateName];
    
    d = sys.d;
    e = sys.e;
    sys.a = blkdiag(sys.a,varargin{i}.a);
    sys.b = blkdiag(sys.b,varargin{i}.b);
    sys.c = blkdiag(sys.c,varargin{i}.c);
    sys.d = blkdiag(d,varargin{i}.d);
    sys.e = blkdiag(e,varargin{i}.e);
    
    sys.y = y;
    sys.u = u;
    sys.StateName = StateName;
    
    if sys.Ts ~= varargin{i}.Ts
        error('Sampling times Ts of models do not match.')
    end
end