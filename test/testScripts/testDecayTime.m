classdef testDecayTime < sssTest
    
    methods (Test)  
        function stableSys(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                warning('off','sss:decayTime:UnstableSys');
                [actTmax]=decayTime(sys);
                if ~isnan(actTmax) % system stable
                    pmin=-log(100)/actTmax; %real(smallest dominant pole)
                    [h]=impulse(ss(pmin,1,1,0),0:actTmax:actTmax); %impulse answer of pmin
                    h=abs(h);
                    expSolution=h(1)*0.01;
                    actSolution=h(2);

                    verification (testCase, actSolution, expSolution);
                    verifyInstanceOf(testCase, actTmax , 'double', 'Instances not matching');
                    verifySize(testCase, actTmax, [1 1], 'Size not matching');
                else
                    warning('off','sss:isstable:EigsFailed');
                    verification(testCase, isstable(sys), 0);
                end
            end
        end
        function unstableSys(testCase)
            warning('off','sss:decayTime:UnstableSys');
            tmax=decayTime(sss(1,1,1)); % unstable system
            verification(testCase,isnan(tmax),true);
        end
        function plausibleTime(testCase)
            %   Confirm through ss/impulse that the decay time computed is
            %   not too far away from the actual one
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                t = decayTime(sys);
                [~,tI] = impulse(ss(sys));
                dt = abs(t-tI(end))/t;
                verifyLessThan(testCase,dt,10) %allow difference to be off by a factor 10
            end
        end
    end
end

function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end