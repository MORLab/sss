function varargout = partition(A,row,col)
% PARTITION - Partitions a matrix into submatrices
%
% Syntax:
%       [A1,A2]             = PARTITION(A,row)
%       [A1,A2]             = PARTITION(A,[],col)
%       [A11,A12,A21,A22]   = PARTITION(A,row,col)
%
% Description:
%       [A1,A2] = PARTITION(A,row) partitions A into two
%       submatrices, the first 1:row and last row+1:end rows, such that 
%       A = [A1;A2].
%
%       [A1,A2] = PARTITION(A,[],col) partitions A into two
%       submatrices, the first 1:col and last col+1:end columns, such that
%       A = [A1,A2].
%
%       [A11,A12,A21,A22] = partition(A,row,col) partitions A into 4 submatrices. 
%       A11 has the dimension (row x col), all other submatrices
%       accordingly, such that A = [A11, A12; A21, A22];
%
% Input Arguments:
%       *Required Input Arguments:*
%       -A:     matrix to be partitioned
%       -row:   number of rows of the first partition
%       *Optional Input Arguments:*
%       -col:   number of columns of the first partition
%
% Output Arguments:
%       -A1,A2:           submatrices resulting from the partition of A
%       -A11,A12,A21,A22: submatrices resulting from the partition of A
%
% See Also:
%       sss
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  06 Feb 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% parse input
if nargin < 3
    %usage partition(A,row)
    col = [];
end

%% partition
if isempty(col)
    A1 = A(1:row,:);
    A2 = A(row+1:end,:);
    
    varargout = {A1,A2};
elseif isempty(row)
    A1 = A(:, 1:col);
    A2 = A(:, col+1:end);
    varargout = {A1,A2};
else
    A11 = A(1:row,1:col);
    A12 = A(1:row,col+1:end);
    A21 = A(row+1:end,1:col);
    A22 = A(row+1:end,col+1:end);
    varargout = {A11,A12,A21,A22};
end
