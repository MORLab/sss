function sys_S = connect(varargin)
% Connects a set of sparse LTI system (sss) by evaluating the names
% of in- and outputs
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% sys_S = connect(varargin);
% Input:        * varargin: sparse state space (sss)-objects
%               * varargin{end}: cell array of output names of the 
%                                closed loop system
%               * varargin{end-1}: cell array of input names of the 
%                                  closed loop system
% Output:       * sys_S: closed loop sparse state space (sss)-object
% ------------------------------------------------------------------
% Authors:      Thomas Emmert (emmert@tfd.mw.tum.de)
% Last Change:  18 Feb 2015
% ------------------------------------------------------------------
% Inspiration: Edwards, J.W., 1976. A FORTRAN program for the analysis of linear continuous and sample-data systems.
% see also: sss/connect_sss sss/append

% Append system to get open loop system
if isa(varargin{end},'cell')
    sys_ap = append(varargin{1:end-2});
    inputname = varargin{end-1};
    outputname = varargin{end};
else
    sys_ap = append(varargin{:});
    inputname = sort(unique(sys_ap.u));
    outputname = sort(unique(sys_ap.y));
end

% Find open loop internal feedbacks of appended system
rows = [];
cols = [];
for k = 1 : length(sys_ap.u)
    col = find(strcmp(sys_ap.y, sys_ap.u(k)));
    
    rows = [rows; k*ones(size(col))];
    cols = [cols; col];
end
vals = ones(size(rows));
K= sparse(rows,cols,vals,sys_ap.m,sys_ap.p);

% Connect internal open loop feedbacks to obtain closed loop model
sys_S = connect_sss(sys_ap, K);

% Truncate internally absorbed and rearrange order of inputs
idx = zeros(length(inputname),1);
for i = 1: length(inputname)
    id = find(strcmp(inputname(i),sys_ap.u));
    if length(id)>1
        warning('Multiple Inputs have the same name.')
    end
    idx(i) = id(1); % TODO: special treatment for one output to multiple inputs (?)
end
u = sys_S.u(idx);
Groups = sys_S.InputGroup;
if not(isempty(Groups))
    for group = fieldnames(Groups)'
        j=1;
        mapping=[];
        for targetidx = 1: length(idx)
            
            sourceIdx = find(sys_S.InputGroup.(char(group))==idx(targetidx));
            if not(isempty(sourceIdx))
                mapping(1,j) = sourceIdx;
                mapping(2,j) = targetidx;
                mapping(3,j) = idx(targetidx); %value
                j= j+1;
            end
        end
        if not(isempty(mapping))
            sys_S.InputGroup.(char(group))(mapping(1,:)) = mapping(2,:);
            sys_S.InputGroup.(char(group)) = sort(unique(sys_S.InputGroup.(char(group))));
        else
            sys_S.InputGroup.(char(group)) = [];
        end
    end
end

sys_S.b = sys_S.b(:,idx);
sys_S.d = sys_S.d(:,idx);
sys_S.u = u;

% Truncate internally absorbed and rearrange order of outputs 
idx = zeros(length(outputname),1);
for i = 1: length(outputname)
    id = find(strcmp(outputname(i),sys_ap.y));
    if length(id)>1
        warning('Multiple Outputs have the same name.')
    end
    idx(i) = id(1);
end
y = sys_S.y(idx);

Groups = sys_S.OutputGroup;
if not(isempty(Groups))
    for group = fieldnames(Groups)'
        j=1;
        mapping=[];
        for targetidx = 1: length(idx)
            
            sourceIdx = find(sys_S.OutputGroup.(char(group))==idx(targetidx));
            if not(isempty(sourceIdx))
                mapping(1,j) = sourceIdx;
                mapping(2,j) = targetidx;
                mapping(3,j) = idx(targetidx); %value
                j= j+1;
            end
        end
        if not(isempty(mapping))
            sys_S.OutputGroup.(char(group))(mapping(1,:)) = mapping(2,:);
            sys_S.OutputGroup.(char(group)) = sort(unique(sys_S.OutputGroup.(char(group))));
        else
            sys_S.OutputGroup.(char(group)) = [];
        end
    end
end

sys_S.c = sys_S.c(idx,:);
sys_S.d = sys_S.d(idx,:);
sys_S.y = y;

end