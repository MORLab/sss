function sys = truncate(sys, idxOut, idxIn)
% TRUNCATE - Truncates a sparse LTI system (sss)
% 
% Syntax:
%       sys = truncate(sys, idxOut, idxIn)
%
% Description:
%       sys = truncate(sys, idxOut, idxIn) truncates the sparse state-space
%       system sys by preserving only the indices defined by idxOut and
%       idxIn
%
% Input Arguments:
%       -sys:    sparse state space (sss)-object
%       -idxOut: indices of outputs to be preserved
%       -idxIn:  indices of inputs to be preserved
%
% Output Arguments:
%       -sys: truncated sparse state-space (sss)-object
%
% Examples:
%       The following code loads the benchmark 'CDplayer' (SSS, MIMO),
%       creates a sss-object and then truncates this system by taking the
%       transfer behaviro from the 2nd input to the 1st output:
%
%> load CDplayer.mat
%> sys=sss(A,B,C); %(SSS, MIMO)
%> TruncatedSys=truncate(sys,1,2); %taking sys12 (SSS, SISO)
%
%       You can compare the results by plotting the frequency response:
%
%> figure; bodemag(sys); %impulse response of the MIMO sss-system
%> figure; bodemag(TruncatedSys); %impulse response of sys12
%
% See Also:
%       sss, plus, minus, mtimes
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
% Authors:      Thomas Emmert
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

idxOut  = string2index(sys.OutputName,idxOut);
idxIn   = string2index(sys.InputName,idxIn);

% truncate matrizes
B = sys.B(:,idxIn);
C = sys.C(idxOut,:);
D = sys.D(idxOut,idxIn);
% truncate IO labels
u = sys.u(idxIn);
y = sys.y(idxOut);
% Assing the values later in order to avoid race condition
sys.B=B; sys.C=C; sys.D=D; sys.u=u; sys.y=y;
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