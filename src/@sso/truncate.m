function sys = truncate(sys, idxOut, idxIn)
% TRUNCATE - Truncates a sso system 
% 
% Syntax:
%       sys = truncate(sys, idxOut, idxIn)
%
% Description:
%       sys = truncate(sys, idxOut, idxIn) truncates the sparse
%       second-order mode sys by preserving only the indices defined by 
%       idxOut and idxIn
%
% Input Arguments:
%       -sys:    sparse second-order (sso)-object
%       -idxOut: indices of outputs to be preserved
%       -idxIn:  indices of inputs to be preserved
%
% Output Arguments:
%       -sys: truncated sso-object
%
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  11 Apr 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

idxOut  = string2index(sys.OutputName,idxOut);
idxIn   = string2index(sys.InputName,idxIn);

% truncate matrizes
B   = sys.B(:,idxIn);
Cp  = sys.Cp(idxOut,:);
Cv  = sys.Cv(idxOut,:);
Df = sys.Df(idxOut,idxIn);

% truncate IO labels
u = sys.u(idxIn);
y = sys.y(idxOut);
% Assing the values later in order to avoid race condition
sys.B=B; sys.Cp=Cp; sys.Cv=Cv; sys.Df=Df; sys.u=u; sys.y=y;
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