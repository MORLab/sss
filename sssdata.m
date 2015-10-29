function [A,B,C,D,E,Ts] = sssdata(sys)
% SSSDATA - Returns all system matrices of an sparse state-space model 
%
% Syntax:
%       [A,B,C,D,E] = SSSDATA(sys);
%       [A,B,C,D,E,Ts] = SSSDATA(sys);
%
% Inputs:
%       -sys: sparse state space (sss)-object
%
%
% Outputs:
%       -A:  system matrix of the sss-object
%       -B:  input matrix of the sss-object
%       -C:  output matrix of the sss-object
%       -D:  feedthrough matrix of the sss-object
%       -E:  descriptor matrix of the sss-object
%       -Ts: sample time  
%
%
% Description:
%       [A,B,C,D,E] = SSSDATA(sys) extracts the system matrices A, B, C, D, E 
%       from the sparse state-space model sys.
%
%       [A,B,C,D,E,Ts] = SSSDATA(sys) also returns the sample time Ts. Other 
%       properties of sys can be accessed using struct-like dot syntax (for
%       example, sys.StateName or sys.x0).
% 
% 
% See also:
%       dssdata, ssdata, sss
%
%
%------------------------------------------------------------------
%   This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order 
%   Reduction and System Analysis Toolbox developed at the Chair of 
%   Automatic Control, Technische Universitaet Muenchen. For updates 
%   and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
%   More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      ??, Maria Cruz Varona 
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  29 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

A = sys.a; B = sys.b; C = sys.C; D = sys.D; E = sys.e;

Ts = sys.Ts;


