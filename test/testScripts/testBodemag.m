classdef testBodemag < sssTest
    
    methods(Test)
        function testCall(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    bodemag(sys,1:100,'r--');
                    bodemag(sys,{10,100});
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
            
            figure; bodemag(sys1,sys2,'r-',sys3,'k--',sys4,sys5);
            figure; bodemag(ss(sys1),ss(sys2),'r-',ss(sys3),'k--',ss(sys4),ss(sys5));
            figure; bodemag(sys1,sys2,'r-',sys3,'k--',{1e-3,1e6},sys4,sys5); % frequency range is passed
        end
        
        function plotFunctionalityMIMO(testCase)
            % verify the correct plot compared to built-in when omega is
            % not passed
            sys1 = sss('CDPlayer');
            sys2 = sss('iss'); sys2 = sys2(1:2,1:2);
            
            figure; bode(sys1,sys2,'r-');
            figure; bode(ss(sys1),ss(sys2),'r-');
            figure; bode(sys1,sys2,'r-',{1e-3,1e6}); % frequency range is passed
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-2,'AbsTol',0.005,...
    'Difference between actual and expected exceeds relative tolerance');
end