classdef testPzmap < sssTest
    % testPzmap - testing of pzmap.m
   
    methods(Test)
        function test1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    % verification of results is done in testZeros and 
                    % testPoles, only test if call works
                    [p,z]=pzmap(sys);
                    pzmap(sys);
                    verifyInstanceOf(testCase, p, 'double');
                    verifyInstanceOf(testCase, z, 'double');
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',0.1e-6,...
    'Difference between actual and expected exceeds relative tolerance');
end