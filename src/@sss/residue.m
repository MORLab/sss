function varargout = residue(varargin)
% RESIDUE - Computes residues, poles and feedthrough of an LTI system
% 
% Syntax:
%       RESIDUE(sys)
%       [r,p,d] = RESIDUE (sys)
%       [r,p,d] = RESIDUE (sys,Opts)
%
% Description:
%       This function computes the pole/residual representation of a
%       rational transfer function (SISO or MIMO). 
%
%       It can either return the residuals themselves, or the residual 
%       directions (low rank factors), useful especially in the MIMO 
%       setting.
%
%       By default, only the first 6 residues corresponding to the
%       eigenvalues with smallest magnitude are computed. This option can
%       be changed by appropriate values in the Opts structure.
%
%       The output r is a cell array of residuals or residual directions
%       depending on the option rType = {'res' (def), 'dir'}.
%
%       For rType = 'res', r is a cell array of dimension nx1 with the
%       residual r{k} for each pole p(k).
%
%       For rType = 'dir' r is a cell array of dimension 1x2 with the
%       output residual matrix Chat = r{1} and the input residual matrix 
%       Bhat = r{2}. The residues can be computed using 
%       r{k} = Chat(:,k)*Bhat(k,:).
%
% Input Arguments:
%       *Required Input Arguments:*
%		-sys: sss-object of the LTI system
%       *Optional Input Arguments:*
%       -Opts:          Structure containing computation options
%           -.rType:    define the output format of the residue r either
%                       as a cell array of dimension nx1 with the residual 
%                       r{k} for each pole p(k) ('res'); or a cell array 
%                       of dimension 1x2 with the output residual matrix 
%                       Chat = r{1} and the input residual matrix Bhat = r{2}
%                       ('dir');
%                       [{'res'} / 'dir']
%			-.eigs:		option to eigs command;
%						[{'SM'} / 'LM' / 'SA' / 'LA' / 'SR' / 'LR' / real or complex scalar]
%           -.nEigs:    number of poles/residues;
%                       [{6} / 'all' / integer]
%
% Output Arguments:
%       -r: cell of residuals (format depends on Opts.rType)
%       -p: eigenvalues
%       -d: feedthrough
%
% Examples:
%		To compute the residuals of a SISO or MIMO system such that 
%       $G(s) = \sum_i(\frac{r_i}{(p_i+s)} + d)$  use
%
%> load building; 
%> sys = sss(A,B,C);
%> Opts.nEig = 'all';
%> [r,p,d] = residue(sys,Opts);
%
%       To get the residue directions instead of the residues, use the
%       option Opts.rType = 'dir'
%
%> load building; 
%> sys = sss(A,B,C);
%> Opts.rType = 'dir';
%> [r,p,d] = residue(sys,Opts);
%
% See Also: 
%		sss, eig, eigs
%
% References:
%		* *[1] Bryson (1994)*, Control of Spacecraft and Aircraft
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
% Authors:      Heiko Panzer, Sylvia Cremer, Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  19 Jan 2016
% Copyright (c) 2015,2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

[varargout{1:nargout}] = sssFunc.residue(varargin{:});

