function z = zero(sys,varargin)
% ZERO - Compute largest invariant zeros of an LTI system
%
% Syntax:
%       z = zero(sys)
%       z = zero(sys,k)
%       z = zero(sys,Opts)
%       z = zero(sys,k,Opts)
%
% Description:
%       z = zero(sys) returns the 6 invariant zeros with largest magnitude
%       of in the column vectors z of the continuous- or discrete-time 
%       dynamic system model sys. The type of the computed zeros can be 
%       specified with the option 'type'.
%
%       z = zero(sys,k) returns the first k zeros of the system.
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
%       Load the benchmark 'CDplayer' (SSS, MIMO) and compute the first 6
%       zeros with largest magnitude:
%
%> load CDplayer.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m))
%> z=zero(sys)
%
% See Also:
%       pzmap, zpk, pole
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
% Last Change:  16 Jun 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% Parse inputs and options
Def.type = 'lm'; %eigs type
Def.sortLm = false; %return first k zeros with 'lm' at sigma

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
    zTemp=cell(sys.p,sys.m);
end

for i=1:sys.p
    for j=1:sys.m
        % call zeros and moments for each siso transfer function

        if strcmp(Opts.type,'lm') || isa(Opts.type,'double')
            
            % use the largest pole instead of 'lm' (usually 'lm' fails)
            if strcmp(Opts.type,'lm')
                try
                    sigma=eigs(sys,1,'lm');
                catch
                    opts.p=4*k;  %double number of lanczos vectors (default: 2*k)
                    sigma=eigs(sys,1,'lm',opts);
                end
            else
                sigma=Opts.type;
            end
            
            % use 1e6 instead of Inf
            if abs(sigma)<1e6
                z=eigs([sys.A,sys.B(:,j);sys.C(i,:),sys.D(i,j)],[sys.E,zeros(sys.n,1);zeros(1,sys.n),0],k,sigma);
                if ~isreal(sigma);
                    % z only contains values with same imaginary sign as temp
                    l=1;
                    while(l<=length(z))
                        if ~isreal(z(l))
                            % add conjugated value if not already in z
                            if min(abs(bsxfun(@minus,[z(1:l-1);z(l+1:end)],conj(z(l)))))>1e-8
                                temp=z((l+1):end);
                                z(l+1)=conj(z(l));
                                if ~isempty(temp)
                                    z(l+2:end+1)=temp;
                                end
                                l=l+1;
                            end
                        end
                        l=l+1;
                    end
                    
                    % remove zeros at infinity
                    z=z(abs(real(z))<1e6);
                    
                    % largest magnitude sorting
                    if strcmp(Opts.type,'lm') || Opts.sortLm
                        tbl=table(-abs(z),z);
                        tbl=sortrows(tbl);
                        z=tbl.z;
                    end
                    z=z(1:k);
                end
            else
                z=eigs([sys.A,sys.B(:,j);sys.C(i,:),sys.D(i,j)],[sys.E,zeros(sys.n,1);zeros(1,sys.n),0],k,-1e6);
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