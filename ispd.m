function [ isPosDef ] = ispd( matrix )
%ISPD  check if the input matrix is positive definite
%
% Syntax:
%   ispd=ISPD(matrix)
%
% Description:
%   This function determines whether the matrix "matrix" is 
%   positive definite or not.
%
% See also:
%   CHOL
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                   -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Jorge Luiz Moreira Silva
% Last Change:  14 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

[~,p] = chol(matrix);
%if and only if the matrix is positive definite, p is zero.
isPosDef=(p==0);
end

