function [A11,A12,A21,A22] = partition(A,row,col)
% PARTITION - Partitions a matrix into four submatrices
%
% Syntax:
%       [A11,A12,A21,A22] = PARTITION(A,row)
%       [A11,A12,A21,A22] = PARTITION(A,row,col)
%
% Description:
%       [A11,A12,A21,A22] = partition(A,row,col) partitions A into 4 submatrices. 
%       A11 has the dimension (row x col), all other submatrices accordingly.
%
%       The input argument |col| is optional. If it is not parsed to the
%       function, then |col = row| is set.
%
% Input Arguments:
%       *Required Input Arguments:*
%           -A:     matrix to be partitioned
%           -row:   number of rows of the submatrix A11
%           -col:   number of columns of the submatrix A11
%
% Output Arguments:
%       -A11,A12,A21,A22: submatrices resulting from the partition of A
%
% Examples:
%       TODO
%
% See Also:
%       TODO
%
% References:
%       TODO
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
% Last Change:  11 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if ~exist('col','var') || isempty(col)
    col = row;
end

A11 = A(1:row,1:col);
if nargout > 1
    A12 = A(1:row,col+1:end);
    A21 = A(row+1:end,1:col);
    A22 = A(row+1:end,col+1:end);
end