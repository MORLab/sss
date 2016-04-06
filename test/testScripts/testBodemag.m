classdef testBodemag < sssTest
    
    methods(Test)
        function testCall(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                bodemag(sys,1:1000,'r--');
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-2,'AbsTol',0.005,...
    'Difference between actual and expected exceeds relative tolerance');
end