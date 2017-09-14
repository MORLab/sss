classdef testBode < sssTest
    
    methods (Test)  
        function mainFunctionality(testCase)
            % verify the correct execution compared to built-in for given
            % frequencies
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                if ~sys_sss.isDae
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
        end
        
        function plotFunctionalitySISO(testCase)
            % verify the correct plot compared to built-in when omega is
            % not passed
            sys1 = sss('building');
            sys2 = sss('beam');
            sys3 = sss('eady');
            sys4 = sss('fom');
            sys5 = sss('iss'); sys5 = sys5(1,1);
            
            figure; bode(sys1,sys2,'r-',sys3,'k--',sys4,sys5);
            figure; bode(ss(sys1),ss(sys2),'r-',ss(sys3),'k--',ss(sys4),ss(sys5));
            figure; bode(sys1,sys2,'r-',sys3,'k--',sys4,sys5,{1e-3,1e6}); % frequency range is passed
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
                if ~sys.isDae
                    w = {10,100};

                    [~,~,omega] = bode(sys,w);
                    verifyEqual(testCase,omega(1),w{1},  'Wrong frequency returned');
                    verifyEqual(testCase,omega(end),w{2},'Wrong frequency returned');
                end
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