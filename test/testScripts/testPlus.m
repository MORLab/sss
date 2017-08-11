classdef testPlus < sssTest
    % testPlus - testing of plus.m
    %
    % Description:
    %   The function plus.m is tested (3 tests) on:
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
    % Authors:      Alessandro Castagnotto, Maria Cruz Varona
    %               Jorge Luiz Moreira Silva
    % Last Change:  05 Nov 2015
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
    
    methods(Test)
        function testplus1(testCase)
            for i=1:length(testCase.sysCell)
                sysSparse=testCase.sysCell{i};
                sys=ss(sysSparse);

                resultSparse = plus(sysSparse, sysSparse);
                result=plus(sys,sys);
                verification(testCase, resultSparse, result);
            end
        end
        function testplus2(testCase)
            initSys2=1;
            for i=1:length(testCase.sysCell)
                if testCase.sysCell{i}.isSiso
                    sys1=testCase.sysCell{i};
                    if initSys2==1
                        sys2=testCase.sysCell{i};
                        initSys2=0;
                    end

                    resultSparse = plus(sys1, sys2);
                    result=plus(ss(sys1),ss(sys2));
                    verification(testCase, resultSparse, result);
                end
            end
        end
        function resultingClass(testCase)
            sys     = testCase.sysCell{1};
            sysRed  = ssRed(sys.A,sys.B,sys.C);
            D       = ones(sys.p,sys.m);
            
            %sparse + ss    -> sss
            sum = sys+ss(sys);
            verifyClass(testCase,sum,'sss');
            
            %ssRed  + sss   -> sss
            sum = sys + sysRed;
            verifyClass(testCase,sum,'sss');
            
            %ssRed  + ssRed -> ssRed
            sum = sysRed + sysRed;
            verifyClass(testCase,sum,'ssRed');
            
            %sss  + D       -> sss
            sum = sys + D;
            verifyClass(testCase,sum,'sss');
            verifyEqual(testCase,sum.D,sparse(sys.D+D)); %correct addition
            verifyEqual(testCase,sys.n,sum.n);   %no dimension increase
            
            %ssRed  + D       -> ssRed
            sum = sysRed + D;
            verifyClass(testCase,sum,'ssRed');
            verifyEqual(testCase,sum.D,sys.D+D); %correct addition
            verifyEqual(testCase,sys.n,sum.n);   %no dimension increase
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