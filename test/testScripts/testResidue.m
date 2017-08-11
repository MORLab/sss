classdef testResidue < sssTest
    % testResidue - testing of residue.m
    %
    % Description:
    %   The function residue.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case just if it runs for the syntax
    %    [r,p,d]=residue(sys) and RESIDUE(sys)
    %
    % ------------------------------------------------------------------
    %   This file is part of sssMOR, a Sparse State Space, Model Order
    %   Reduction and System Analysis Toolbox developed at the Institute
    %   of Automatic Control, Technische Universitaet Muenchen.
    %   For updates and further information please visit www.rt.mw.tum.de
    %   For any suggestions, submission and/or bug reports, mail us at
    %                     -> morlab@rt.mw.tum.de <-
    % ------------------------------------------------------------------
    % Authors:      Alessandro Castagnotto
    %               Jorge Luiz Moreira Silva
    % Last Change:  26 Out 2015
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
   
    methods(Test)
        function fullResidue(testCase)
            Opts.nEigs = 'all';
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                
                [r,p,d]=residue(sys,Opts);
                l = eig(sys);
                
                verifyEqual(testCase,cplxpair(p),cplxpair(l),'RelTol',1e-3);
                verifyEqual(testCase,d,sys.d,'RelTol',1e-3);
                verifySize(testCase,r,[sys.n,1])
            end
        end
        
        function sparseResidueDefault(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                
                [r,p,d]=residue(sys);
                l = eigs(sys,6,'sm');
                
                verifyEqual(testCase,cplxpair(p),cplxpair(l),'RelTol',1e-3);
                verifyEqual(testCase,d,sys.d,'RelTol',1e-3);
                verifySize(testCase,r,[6,1])
            end
        end
        function sparseResidueCustom(testCase)
            Opts.eigs   = 'lm';
            Opts.nEigs  = 10;
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                
                [r,p,d]=residue(sys,Opts);
                l = eigs(sys,Opts.nEigs,Opts.eigs);
                
                verifyEqual(testCase,cplxpair(p),cplxpair(l),'RelTol',1e-3);
                verifyEqual(testCase,d,sys.d,'RelTol',1e-3);
                verifySize(testCase,r,[Opts.nEigs,1])
            end
        end
    end
end
