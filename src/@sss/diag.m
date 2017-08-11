function varargout = diag(varargin)
% DIAG - Transforms an LTI system to (block)-diagonal representation
%
% Syntax:
%       sysd = DIAG(sys)
%
% Description:
%       sysd = DIAG(sys) transforms the sparse state-space model sys to a 
%       (block)diagonal representation. The resulting diagonal state-space
%       representation sysd is always of type sss (E = I), both in case
%       where the original model sys is sss (E = I) or dss (E ~= I, E invertible). 
%
%       If sys has real eigenvalues, then the diagonal matrix sysd.A
%       contains the real eigenvalues in the diagonal. If sys has complex
%       conjugate eigenvalues with real part delta and imaginary part omega, 
%       then the diagonal matrix sysd.A has a block-diagonal structure:
%
%       |sysd.A = [delta1  omega1    0      0;
%                 -omega1 delta1    0      0;
%                    0       0   delta2  omega2;
%                    0       0   -omega2 delta2]|.
%       
%       During the diagonalization, the C-vector is normalized to contain
%       ones.
%
%       //Note: this function performs dense computations and is hence
%       suitable only for mid-sized problems
%
% Input Arguments:
%       -sys: an sss-object containing the LTI system
%
% Output Arguments:
%       -sysd: system in diagonal form
%
% Examples:
%       To compute the diagonal state-space realization of the benchmark
%       "building" (SSS, SISO) use:
%
%> load building.mat
%> sys = sss(A,B,C)
%> sysd = diag(sys)
%> [Ad,Bd,Cd] = ssdata(sysd);
%
%       You can compare the sparsity patterns of |A| and |Ad|:
%
%> figure;
%> subplot(1,2,1); spy(sys.A); title('Sparsity pattern of A');
%> subplot(1,2,2); spy(Ad); title('Sparsity pattern of Ad');
%
%       DIAG also supports SIMO, MISO and MIMO as well as DSSS systems.
%       The following code loads the benchmark 'rail_1357' (DSSS, MIMO)
%       and transforms it into a diagonal representation with E=I (SSS,
%       MIMO):
% 
%> load rail_1357.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m),E)
%> sysd = diag(sys)
%> [Ad,Bd,Cd,Dd] = ssdata(sysd);
%
% See Also: 
%        eig, residue
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
% Authors:      Heiko Panzer, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  04 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

[varargout{1:nargout}] = sssFunc.diag(varargin{:});
