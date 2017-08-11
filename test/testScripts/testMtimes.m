classdef testMtimes < sssTest
    % testMtimes - testing of mtimes.m
%
% Description:
%   The function rk.m is tested (3 tests) on:
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
%               Jorge Luiz Moreira Silva
% Last Change:  26 Out 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

    methods(Test)
        function testMTimes1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
            
                resultSparse = mtimes(sys, sys');
                result=mtimes(ss(sys),ss(sys'));
                verification(testCase, resultSparse, result);
            end
        end
        function testMTimes2(testCase)
            initSys2=1;
            for i=1:length(testCase.sysCell)
                if testCase.sysCell{i}.isSiso
                    sys1=testCase.sysCell{i};
                    if initSys2==1
                        sys2=testCase.sysCell{i};
                        initSys2=0;
                    end

                    resultSparse = mtimes(sys1, sys2);
                    result=mtimes(ss(sys1),ss(sys2));
                    verification(testCase, resultSparse, result);
                    sys2=sys1;
                end
            end
        end
        function resultingClass(testCase)
            sys     = testCase.sysCell{1};
            sysRed  = ssRed(sys.A,sys.B,sys.C);
            
            %sparse * ss    -> sss
            prod = sys*ss(sys);
            verifyClass(testCase,prod,'sss');
            
            %ssRed  * sss   -> sss
            prod = sys * sysRed;
            verifyClass(testCase,prod,'sss');
            
            %ssRed  * ssRed -> ssRed
            prod = sysRed * sysRed;
            verifyClass(testCase,prod,'ssRed');
            
        end
    end
end

function [] = verification(testCase, actSolution, expSolution, m)
verifyEqual(testCase, full(actSolution.A),  full(expSolution.A),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
verifyEqual(testCase,  full(actSolution.B),  full(expSolution.B),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
verifyEqual(testCase,  full(actSolution.C),  full(expSolution.C),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
verifyEqual(testCase,  full(actSolution.D),  full(expSolution.D),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
end