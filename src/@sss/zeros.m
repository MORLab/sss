function z = zeros(sys,varargin)
% ZEROS - Compute largest invariant zeros of an LTI system
%
% Syntax:
%       z = zeros(sys)
%       z = zeros(sys,k)
%       z = zeros(sys,Opts)
%       z = zeros(sys,k,Opts)
%
% Description:
%       z = zeros(sys) returns the 6 invariant zeros with largest magnitude
%       of in the column vectors z of the continuous- or discrete-time 
%       dynamic system model sys. The type of the computed zeros can be 
%       specified with the option 'type'.
%
%       z = zeros(sys,k) returns the first k zeros of the system.
%
%//Note: The calculation of the invariant zeros is only defined for systems
%       with the same number of inputs and outputs (m=p). That means that if
%       zpk is called with a system with m~=p, then z = [ ].
%
% Input Arguments:
%       -sys:      an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -k:     number of computed zeros
%       -Opts:  structure with execution parameters
%			-.type:  eigs type;
%						[{'lm'} / 'sm' / 'la' / 'sa']
%
% Output Arguments:
%       -z: vector containing invariant zeros
%
% Examples:
%       Load the benchmark 'rail_1357' (DSSS, MIMO) and compute the first 6
%       zeros with largest magnitude:
%
%> load rail_1357.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m),E)
%> z=zeros(sys);
%
% See Also:
%       pzmap, zpk, poles
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
% Authors:      Alessandro Castagnotto, Maria Cruz Varona,
%               Lisa Jeschek
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
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

if sys.isDae
    error('Zeros does not work with DAE systems yet.');
elseif ~sys.isSiso
    error('Zeros works only with SISO systems yet.');
end

% zero function
if sys.m==sys.p    
    if ~exist('k','var')
        k=6;
    end
    
    if strcmp(Opts.type,'sm')
        Opts.type=1e-2;
        E22=0;
    elseif strcmp(Opts.type,'lm')
        E22=1e-16;
    else %sigma
        E22=0;
    end
    
    z=eigs([sys.A,sys.B;sys.C,sys.D],[sys.E,zeros(sys.n,sys.m);zeros(sys.p,sys.n),E22],k,Opts.type);

    % remove zeros at infinity
    z=z(abs(real(z))<1e6);
    
    % compare results of eig and eigs
    zEig=eig(full([sys.A,sys.B;sys.C,sys.D]),[full(sys.E),zeros(sys.n,sys.m);zeros(sys.p,sys.n),zeros(sys.p,sys.m)]);
    
    % remove zeros at infinity
    zEig=zEig(abs(real(zEig))<1e6);
    
    if strcmp(Opts.type,'lm')
        tbl=table(-abs(z),z);
        tbl=sortrows(tbl);
        z=tbl.z;
        tbl=table(-abs(zEig),zEig);
        tbl=sortrows(tbl);
        zEig=tbl.zEig;
    else
        tbl=table(abs(z),z);
        tbl=sortrows(tbl);
        z=tbl.z;
        tbl=table(abs(zEig),zEig);
        tbl=sortrows(tbl);
        zEig=tbl.zEig;
    end
    
    % remove single complex element (for comparison)
    if abs(imag(sum(z)))>1e-12
       z(abs(imag(z)-imag(sum(z)))<1e-16)=[];
    end

    % display all values if difference is big
    if norm(zEig(1:size(z,1))-z)>1
        disp([z,zEig(1:size(z,1))]);
    end
    
    % display difference
    disp('norm(z_eig-z_eigs):');
    disp(norm(zEig(1:size(z,1))-z))
    
else
    z=zeros(0,1);
end

end