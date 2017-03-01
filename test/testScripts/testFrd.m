classdef testFrd < sssTest
    
    methods (Test)
        function generalFunctionality(testCase)
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                if ~sys_sss.isDae
                    sys_ss=ss(sys_sss);

                    [expSolution, ~, omega]=bode(sys_ss);

                    frdobj=frd(sys_sss,omega);
                    actSolution=abs(frdobj.responseData);

                    verification(testCase, actSolution, expSolution);
                    verifyInstanceOf(testCase, frdobj, 'frd', 'Instances not matching');

                    frdobj=frd(sys_sss);
                    verifyInstanceOf(testCase, frdobj, 'frd', 'Instances not matching');
                end
            end
        end
    end
end

function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end