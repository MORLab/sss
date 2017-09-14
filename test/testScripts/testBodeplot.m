classdef testBodeplot < sssTest
    
    methods(Test)
        function testBodeplot1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    h=bodeplot(sys,1:100,'r--');
                    setoptions(h,'PhaseVisible','off');
                    verifyInstanceOf(testCase, h, 'handle','Instances not matching');
                    bodeplot(sys,{10,100});
                end
            end
        end
        
        function plotFunctionalitySISO(testCase)
            % verify the correct plot compared to built-in when omega is
            % not passed
            sys1 = sss('building');
            sys2 = sss('beam');
            sys3 = sss('eady');
            sys4 = sss('fom');
            sys5 = sss('iss'); sys5 = sys5(1,1);
            
            figure; bodeplot(sys1,sys2,'r-',sys3,'k--',sys4,sys5);
            figure; bodeplot(ss(sys1),ss(sys2),'r-',ss(sys3),'k--',ss(sys4),ss(sys5));
            figure; bodeplot(sys1,sys2,'r-',sys3,'k--',sys4,sys5,{1e-3,1e6}); % frequency range is passed
        end
        
        function plotFunctionalityMIMO(testCase)
            % verify the correct plot compared to built-in when omega is
            % not passed
            sys1 = sss('CDPlayer');
            sys2 = sss('iss'); sys2 = sys2(1:2,1:2);
            
            figure; bodeplot(sys1,sys2,'r-');
            figure; bodeplot(ss(sys1),ss(sys2),'r-');
            figure; bodeplot(sys1,sys2,'r-',{1e-3,1e6}); % frequency range is passed
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-2,'AbsTol',0.005,...
    'Difference between actual and expected exceeds relative tolerance');
end