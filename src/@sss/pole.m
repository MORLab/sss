function p = pole(sys,varargin)
% POLE - Compute largest poles of an LTI system
%
% Syntax:
%       p = pole(sys)
%       p = pole(sys,k)
%       p = pole(sys,Opts)
%       p = pole(sys,k,Opts)
%
% Description:
%       p = pole(sys) returns the 6 poles with largest magnitude of in the
%       column vectors z of the continuous- or discrete-time dynamic system
%       model sys. The type of the computed poles can be specified with the
%       option 'type'.
%
%       p = pole(sys,k) returns the first k poles of the system.
%
% Input Arguments:
%       -sys:      an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -k:     number of computed poles
%       -Opts:  structure with execution parameters
%			-.type:  eigs type;
%						[{'lm'} / 'sm' / 'la' / 'sa']
%
% Output Arguments:
%       -p: vector containing the system poles
%
% Examples:
%       Load the benchmark 'rail_1357' (DSSS, MIMO) and compute the first 6
%       poles with largest magnitude:
%
%> load rail_1357.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m),E)
%> p=pole(sys)
%
% See Also:
%       pzmap, zpk, zero
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
% Authors:      Alessandro Castagnotto, Maria Cruz Varona,
%               Lisa Jeschek
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  19 Apr 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% Parse inputs and options
Def.type = 'lm'; %eigs type

for i=1:length(varargin)
    if isa(varargin{i},'double')
        k=varargin{i};
    elseif isa(varargin{i},'struct')
        Opts=varargin{i};
    end
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

if ~exist('k','var')
    k=6;
end

if ~sys.isDae
    try
        [~,p,flag] = eigs(sys.A,sys.E,k,Opts.type);
        if flag~=0
            error('Not all eigenvalues converged');
        end
        p=diag(p);
    catch
        opts.p=4*k; %double number of lanczos vectors (default: 2*k)
        p = eigs(sys.A,sys.E,k,Opts.type,opts);
    end
else
    error('Poles does not work with DAE systems yet.');
end

% ensure column vector
if size(p,1)<size(p,2)
    p=transpose(p);
end

end