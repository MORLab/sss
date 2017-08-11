classdef testMinus < sssTest
    % testMinus - testing of minus.m
    %
    % Description:
    %   The function minus.m is tested (3 tests) on:
    %    + combination of two benchmark-systems.
    %    + combination of two random-systems that are equal.
    %    + combination of two random-systems that are different.
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
    % Last Change:  13 Apr 2017
    % Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
    
    methods(Test)
        function benchmarkTest(testCase)
            for i=1:length(testCase.sysCell)
                sysSparse=testCase.sysCell{i};
                sys=ss(sysSparse);

                resultSparse = minus(sysSparse, sysSparse);
                result=minus(sys,sys);
                verification(testCase, resultSparse, result);
            end
        end
        function benchmarkTest2(testCase)
            initSys2=1;
            for i=1:length(testCase.sysCell)
                if testCase.sysCell{i}.isSiso
                    sys1=testCase.sysCell{i};
                    if initSys2==1
                        sys2=testCase.sysCell{i};
                        initSys2=0;
                    end

                    resultSparse = minus(sys1, sys2);
                    result = minus(ss(sys1),ss(sys2));
                    verification(testCase, resultSparse, result);
                end
            end
        end
        function resultingClass(testCase)
            sys     = testCase.sysCell{1};
            sysRed  = ssRed(sys.A,sys.B,sys.C);
            D       = ones(sys.p,sys.m);
            
            %sparse - ss    -> sss
            diff = sys-ss(sys);
            verifyClass(testCase,diff,'sss');
            
            %ssRed  - sss   -> sss
            diff = sys - sysRed;
            verifyClass(testCase,diff,'sss');
            
            %ssRed  - ssRed -> ssRed
            diff = sysRed - sysRed;
            verifyClass(testCase,diff,'ssRed');
            
            %sss  - D       -> sss
            diff = sys - D;
            verifyClass(testCase,diff,'sss');
            verifyEqual(testCase,diff.D,sparse(sys.D-D)); %correct addition
            verifyEqual(testCase,sys.n,diff.n);   %no dimension increase
            
            %ssRed  - D       -> ssRed
            diff = sysRed - D;
            verifyClass(testCase,diff,'ssRed');
            verifyEqual(testCase,diff.D,sys.D-D); %correct addition
            verifyEqual(testCase,sys.n,diff.n);   %no dimension increase
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, full(actSolution.A),  full(expSolution.A),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
verifyEqual(testCase,  full(actSolution.B),  full(expSolution.B),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
verifyEqual(testCase,  full(actSolution.C),  full(expSolution.C),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
verifyEqual(testCase,  full(actSolution.D),  full(expSolution.D),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
end