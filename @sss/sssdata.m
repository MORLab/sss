function [A,B,C,D,E] = sssdata(sys_sss)
% Returns all system matrices in one expression
% ------------------------------------------------------------------
% [A,B,C,D,E] = sssdata(sys_sss);
% Input:        * sys_sss: sparse state space (sss)-object
% Output:       * A: system matrix of the sss-object
%               * B: input matrix of the sss-object
%               * C: output matrix of the sss-object
%               * D: feedthrough matrix of the sss-object
%               * E: descriptor matrix
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      ??
% Last Change:  22 Sep 2015
% ------------------------------------------------------------------

A = sys_sss.a; B = sys_sss.b; C = sys_sss.C; D = sys_sss.D; E = sys_sss.e;

