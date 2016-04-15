classdef testBode < sssTest
    
    methods (Test)  
        function mainFunctionality(testCase)
            % verify the correct execution compared to built-in for given
            % frequencies
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                t=1:100:1000;
                
                [actMag, actPhase, actOmega]=bode(sys_sss,t);
                [expMag, expPhase, expOmega]=bode(ss(sys_sss),t);
                
                %Phase between 0? to 360?
                for j=1:length(actPhase)
                    %if actPhase(:,:,j)<0
                    actPhase(:,:,j)=actPhase(:,:,j)-floor(actPhase(:,:,j)/360)*360;
                    %end
                    %if expPhase(:,:,j)<0
                    expPhase(:,:,j)=expPhase(:,:,j)-floor(expPhase(:,:,j)/360)*360;
                    %end
                end
                
                actSolution={actMag, actPhase, actOmega};
                expSolution={expMag, expPhase, expOmega};
                
                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actMag , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actPhase , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actOmega , 'double', 'Instances not matching');
                verifySize(testCase, actMag, size(expMag), 'Size not matching');
                verifySize(testCase, actPhase, size(expPhase), 'Size not matching');
                verifySize(testCase, actOmega, size(expOmega), 'Size not matching');
            end
        end
        function outputFunctionality(testCase)
            % FRD-Object functionality
            % 1) return an error if more than one system was passed and
            %    output variables are defined
            % 2) return a magnitude array if no option is passed
            
            for i=1:length(testCase.sysCell)
                sys  = testCase.sysCell{i};
 
                w = 1:100:1000;
                
                verifyError(testCase, @() triggerError(sys,sys),...
                    'sss:bode:RequiresSingleModelWithOutputArgs');
                
                mag = bode(sys,w);
                    verifyClass(testCase,mag,'double');
                    
            end
        end
        function inputFunctionality(testCase)
            % input of frequnecy range
            
            for i=1:length(testCase.sysCell)
                sys  = testCase.sysCell{i};
                w = {10,100};

                [~,~,omega] = bode(sys,w);
                verifyEqual(testCase,omega(1),w{1},  'Wrong frequency returned');
                verifyEqual(testCase,omega(end),w{2},'Wrong frequency returned');
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end

function triggerError(varargin)
     bla = bode(varargin{:});
end