function varargout = zeros(varargin)
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
%> load CDplayer.mat
%> p = size(C,1); m = size(B,2);
%> sys = sss(A,B,C,zeros(p,m))
%> z=zeros(sys)
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
% Last Change:  16 Jun 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

[varargout{1:nargout}] = sss.zeros(varargin{:});