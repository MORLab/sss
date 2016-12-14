function zpkData = zpk(sys,varargin)
% ZPK - Compute largest poles and zeros or zpk object of an LTI system
%
% Syntax:
%       zpkData = ZPK(sys)
%       zpkData = ZPK(sys,kP,typeP,kZ,typeZ)
%       zpkData = ZPK(sys,k,typeP,typeZ)
%       zpkData = ZPK(sys,kP,kZ,type)
%       zpkData = ZPK(sys,kP,kZ)
%       zpkData = ZPK(sys,typeP,typeZ)
%       zpkData = ZPK(sys,k,type)
%       zpkData = ZPK(sys,k)
%       zpkData = ZPK(sys,type)
%
% Description:
%       [p,z] = zpk(sys) returns the 6 system poles and invariant zeros 
%       with largest magnitude of in the column vectors p and z of the 
%       continuous- or discrete-time dynamic system model sys. The type of 
%       the computed poles and zeros can be specified with the options 
%       'typeP' and 'typeZ'.
%
%       [p,z] = zpk(sys,k) returns the first k poles and zeros of the system.
%       
%       If the option 'zpk' is true, a zpk-object is returned instead of
%       the poles and zeros.
%
%//Note: If the system is MIMO, the zeros are computed for all combination
%       of inputs and outputs and z is returned in a cell array.
%
% Input Arguments:
%       -sys:      an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -k:     number of computed poles and zeros
%       -Opts:  structure with execution parameters
%			-.zpk:  return zpk object;
%						[{0} / 1]
%			-.typeP: eigs type of poles
%						[{'lm'} / 'sm' / 'la' / 'sa']
%			-.typeZ: eigs type of zeros
%						[{'lm'} / 'sm' / 'la' / 'sa']
%
% Output Arguments:
%       -p: vector containing poles 
%       -z: vector/cell array containing invariant zeros
%
% Examples:
%       Create a random descriptor model (DSSS, SISO) and compute the poles
%       and zeros.
%
%> A = randn(500,500); B = randn(500,1); C = randn(1,500); D = zeros(1,1);
%> E = randn(500,500);
%> sys = dss(A,B,C,D,E);
%> sysSss = sss(sys);
%> [p,z]=zpk(sysSss)
%
%       Load the benchmark 'CDplayer' (SSS, MIMO) and use zpk:
%
%> load CDplayer.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m))
%> [p,z]=zpk(sys)   %computing poles and zeros
%> zpkData=zpk(sys) %computing zpkData (numerator and denominator data)
%
% See Also:
%       pzmap, zeros, poles
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
% Last Change:  16 Jun 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% Parse inputs
Opts.typeZ='lm';Opts.typeP='lm';kZ = 6;kP = 6;   % default values
try
    if nargin == 2
        if isInt(varargin{1})
            kP = varargin{1};
            kZ = varargin{1};
        elseif isOpts(varargin{1})
            Opts.typeP = varargin{1};
            Opts.typeZ = varargin{1};
        else error(''); end
    elseif nargin == 3
        if isOpts(varargin{1}) && isOpts(varargin{2})
            Opts.typeP = varargin{1};
            Opts.typeZ = varargin{2};
        elseif isInt(varargin{1}) && isInt(varargin{2})
            kP = varargin{1};
            kZ = varargin{2};
        elseif isInt(varargin{1}) && isOpts(varargin{2})
            kP = varargin{1};
            kZ = varargin{1};
            Opts.typeP = varargin{2};
            Opts.typeZ = varargin{2};
        else error(' '); end
    elseif nargin == 4
        if isInt(varargin{1}) && isOpts(varargin{3})
            kP = varargin{1};
            if isInt(varargin{2})
                kZ = varargin{2};
                Opts.typeZ = varargin{3};
                Opts.typeP = varargin{3};
            elseif isOpts(varargin{2})
                kZ = kP;
                Opts.typeP = varargin{2};
                Opts.typeZ = varargin{3};
            else error(''); end
        else error(' '); end
    elseif nargin == 5
        if isInt(varargin{1}) && isOpts(varargin{2}) && ... 
           isInt(varargin{3}) && isOpts(varargin{4})
            kP = varargin{1};
            kZ = varargin{2};
            Opts.typeP = varargin{3};
            Opts.typeZ = varargin{4};
        else error(' '); end
    end
catch ex
   error('Arguments have the wrong format. Type "help sss/zpk" for more information.')
end

%% Calculate poles and zeros
pTemp=poles(sys,kP,struct('type',Opts.typeP));

p=cell(sys.m,sys.p);
z=cell(sys.m,sys.p);
c=zeros(sys.m,sys.p);

if strcmp(Opts.typeZ,'lm')
    Opts.type=max(pTemp);
    Opts.sortLm=true;
end

% remove single complex element
if ~isreal(pTemp)
    temp=pTemp(abs(imag(pTemp)-imag(sum(pTemp)))<1e-16);
    if ~isempty(temp) && any(abs(imag(temp))>1e-16)
        pTemp(end+1)=temp;
    end
end

for i=1:sys.m
    for j=1:sys.p
        % call zeros and moments for each siso transfer function
        tempSys=sss(sys.A,sys.B(:,j),sys.C(i,:),sys.D(i,j),sys.E);
        zTemp=zeros(tempSys,kZ,Opts);

        % remove not converged eigenvalues
        pTemp(isnan(pTemp))=[];
        zTemp(isnan(zTemp))=[];

        % avoid infinity
        pTemp(abs(pTemp)>1e6)=1e6*sign(real(pTemp(abs(pTemp)>1e6)));
        zTemp(abs(zTemp)>1e6)=1e6*sign(real(zTemp(abs(zTemp)>1e6)));

        % add second element of single complex element
        if ~isreal(zTemp)
            temp=zTemp(abs(imag(zTemp)-imag(sum(zTemp)))<1e-16);
            if ~isempty(temp)
                for k=1:size(temp,1)
                   zTemp(end+1)=conj(temp(k));
                end
            end
        end

        p{i,j}=pTemp;
        z{i,j}=zTemp;

        % gain c is the first nonzero markov parameter
        ctemp=moments(tempSys,Inf,2);
        c(i,j)=ctemp(:,:,2);
    end
end
 
zpkData=zpk(z,p,c);
zpkData.Name=sys.Name;
    
    
function tf = isInt(A)
    if isnumeric(A) && mod(A,1) == 0
        tf = 1;
    else
        tf = 0;
    end

function tf = isOpts(A)
    if ischar(A) && ismember(A,{'lm','sm','la','sa'})
       tf = 1;
    else
       tf = 0;
    end