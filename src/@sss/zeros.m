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
end

if ~exist('k','var')
    k=6;
end

if strcmp(Opts.type,'sm')
    Opts.type=1e-2;
end

% zero function
if ~sys.isSiso
    zTemp=cell(sys.m,sys.p);
end

for i=1:sys.m
    for j=1:sys.p
        % call zeros and moments for each siso transfer function

        if strcmp(Opts.type,'lm')
            try
                % sigma & E22=0, does not work with all systems,
                % result does not contain wrong values, but some values may
                % not be included
                temp=eigs(sys,1,'lm'); %sigma = max magnitude pole
                if abs(temp)<1e6 % not infinity
                    z=eigs([sys.A,sys.B(:,j);sys.C(i,:),sys.D(i,j)],[sys.E,zeros(sys.n,1);zeros(1,sys.n),0],k,temp);
                    if ~isreal(temp);
                        % finds all values with the same sign of the
                        % imaginary part like temp, add conjugated values to
                        % result vector
                        z=[z;conj(z)];
                    end
                else % infinity
                    % eigs fails with Inf or very big values, try infinity
                    % threshold of 1e6 instead
                    z=eigs([sys.A,sys.B(:,j);sys.C(i,:),sys.D(i,j)],[sys.E,zeros(sys.n,1);zeros(1,sys.n),0],k,-1e6);
                end
            catch
                % 1e-16 instead of E22, works for all not-dae system, but
                % result may contain wrong values
                z=eigs([sys.A,sys.B(:,j);sys.C(i,:),sys.D(i,j)],[sys.E,zeros(sys.n,1);zeros(1,sys.n),1e-16],k,Opts.type);
            end
        else
            z=eigs([sys.A,sys.B(:,j);sys.C(i,:),sys.D(i,j)],[sys.E,zeros(sys.n,1);zeros(1,sys.n),0],k,Opts.type);
        end

        % remove zeros at infinity
        z=z(abs(real(z))<1e6); 
        if ~sys.isSiso
            zTemp{i,j}=z;
        end
    end
end
if ~sys.isSiso
    z=zTemp;
end

end