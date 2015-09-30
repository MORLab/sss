function sys = append(varargin)
% Appends a set of sparse LTI system (sss) 
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% sys = append(varargin)
% Input:        * varargin: sparse state space (sss)-objects
% Output:       * sys_S: appended (open loop) sparse state space 
%                        (sss)-object
% ------------------------------------------------------------------
% Authors:      Thomas Emmert (emmert@tfd.mw.tum.de)
% Last Change:  20 Feb 2015
% ------------------------------------------------------------------
% see also sss/connect

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