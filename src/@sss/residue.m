function [r,p,d] = residue(sys, Opts)
% RESIDUE - Computes residues, poles and feedthrough of an LTI system
% 
% Syntax:
%       RESIDUE(sys)
%       [r] = RESIDUE(sys)
%       [r,p] = RESIDUE (sys)
%       [r,p,d] = RESIDUE(sys,Opts)
%
% Description:
%       This function computes the pole/residual representation of a
%       rational transfer function (SISO or MIMO). 
%
%       It can either return the residuals themselves, or the residual 
%       directions (low rank factors), useful especially in the MIMO 
%       setting.
%
%       The computation requires the complete eigendecomposition of the
%       system, so for large scale problems it might take a while.
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
%> [r,p,d] = residue(sys);
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
%		sss, eig
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
%                   -> <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Sylvia Cremer, Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  28 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

Def.rType = 'res';

if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end   

% are poles and residues already available?
if ~isempty(sys.poles) && ~isempty(sys.residues)
    p=sys.poles;
    r=sys.residues;
    d=sys.D;
    return
end

%perform eigen-decomposition of system
try
    [T,J] = eig(sys);
catch err
    error('Computation of the eigenvalues and eigenvectors failed with message:%s',err.message);
end

% transform system to diagonal form
p=diag(J).';
if issparse(T)
    rcondNumber = 1/condest(T);
else
    rcondNumber=rcond(T);
end
if rcondNumber<eps
    warning('Matrix of eigenvectors is close to singular or badly scaled. Results may be inaccurate. RCOND =',num2str(rcondNumber));
end
B=(sys.E*T)\sys.B;
C=sys.C*T;
d=sys.D;

% calculate residues
r = cell(1,sys.n);
for i=1:sys.n
    r{i} = full(C(:,i)*B(i,:));
end

% store results to caller workspace
sys.residues=r;
sys.poles=p;
if inputname(1)
    assignin('caller', inputname(1), sys);
end

% return the residual directions instead of the residuals
if strcmp(Opts.rType,'dir')
    clear r  
    r = {C, B};  
end
return
