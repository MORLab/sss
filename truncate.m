function sys = truncate(sys, idxOut, idxIn)
% truncates a sparse LTI system (sss)
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% sys = truncate(sys, idxOut, idxIn);
% Input:        * sys:    sparse state space (sss)-object
%               * idxOut: indices of outputs to be preserved
%               * idxIn:  indices of inputs to be preserved
% Output:       * sys_S:  appended (open loop) sparse state space 
%                         (sss)-object
% ------------------------------------------------------------------
% Authors:      Thomas Emmert (emmert@tfd.mw.tum.de)
% Last Change:  18 Apr 2015
% ------------------------------------------------------------------

idxOut = string2index(sys.OutputName,idxOut);
idxIn = string2index(sys.InputName,idxIn);

% truncate matrizes
sys.B = sys.B(:,idxIn);
sys.C = sys.C(idxOut,:);
sys.D = sys.D(idxOut,idxIn);
% truncate IO labels
% sys.u = sys.u(idxIn);
% sys.y = sys.y(idxOut);
% truncate the groups:
InGroups = sys.InputGroup;
if not(isempty(InGroups))
    for Groupc = fieldnames(InGroups)'
        if not(isempty(Groupc))
            Group= char(Groupc);
            sys.InputGroup.(Group) = find(ismember(idxIn, sys.InputGroup.(Group)));
        end
    end
end
OutGroups = sys.OutputGroup;
if not(isempty(OutGroups))
    for Groupc = fieldnames(OutGroups)'
        if not(isempty(Groupc))
            Group= char(Groupc);
            sys.OutputGroup.(Group) = find(ismember(idxOut, sys.OutputGroup.(Group)));
        end
    end
end

end

function idx = string2index(signalNames,idx)
if ischar(idx)
    idx = {idx};
end
if iscell(idx)    
    index = [];
    for i = 1:length(idx)
        if strcmp(idx{i},':')
            idx{i} = '.*';
        else
            idx{i} = ['^' idx{i} '$'];
        end
        index = [index find(cellfun(@(x) ~isempty(regexp(x,idx{i})),signalNames))'];
    end
    idx = unique(index,'stable');
end

end