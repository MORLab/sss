function [varargout] = zpk(varargin)
% ZPK - Compute largest poles and zeros or zpk object of an LTI system
%
% Syntax:
%       [p,z] = ZPK(sys)
%       [p,z] = ZPK(sys,k)
%       zpkData = ZPK(sys,Opts)
%       zpkData = ZPK(sys,k,Opts)
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

[varargout{1:nargout}] = sss.zpk(varargin{:});