classdef testBodeplot < sssTest
    
    methods(Test)
        function testBodeplot1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                h=bodeplot(sys,1:100,'r--');
                setoptions(h,'PhaseVisible','off');
                verifyInstanceOf(testCase, h, 'handle','Instances not matching');
                bodeplot(sys,{10,100});
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-2,'AbsTol',0.005,...
    'Difference between actual and expected exceeds relative tolerance');
end