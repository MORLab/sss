classdef testLyapchol < sssTest
    methods(Test)
        function testLyapchol1(testCase) 
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    [actR,actL]=lyapchol(sys);

                    if sys.isDescriptor
                        actResR=norm(sys.A*(actR'*actR)'*sys.E'+sys.E*(actR'*actR)'*sys.A'+sys.B*sys.B');
                        actResL=norm(sys.A'*(actL'*actL)'*sys.E+sys.E'*(actL'*actL)'*sys.A+sys.C'*sys.C);
                    else
                        actResR=(sys.A*(actR'*actR)+(actR'*actR)*sys.A'+sys.B*sys.B');
                        actResL=(sys.A'*(actL'*actL)+(actL'*actL)*sys.A+sys.C'*sys.C);
                    end

                    actSolution={actResR, actResL};
                    expSolution={zeros(size(actResR)),zeros(size(actResL))};
                    verification(testCase, actSolution, expSolution);
                end
            end
        end
        function testLyapcholOpts(testCase)
            Opts.method='adi';
            Opts.q=120;
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && sys.n>300
                    [actR,actL]=lyapchol(sys,Opts);

                    if sys.isDescriptor
                        actResR=norm(sys.A*(actR'*actR)'*sys.E'+sys.E*(actR'*actR)'*sys.A'+sys.B*sys.B');
                        actResL=norm(sys.A'*(actL'*actL)'*sys.E+sys.E'*(actL'*actL)'*sys.A+sys.C'*sys.C);
                    else
                        actResR=(sys.A*(actR'*actR)+(actR'*actR)*sys.A'+sys.B*sys.B');
                        actResL=(sys.A'*(actL'*actL)+(actL'*actL)*sys.A+sys.C'*sys.C);
                    end

                    actSolution={actResR, actResL};
                    expSolution={zeros(size(actResR)),zeros(size(actResL))};
                    verification(testCase, actSolution, expSolution);
                end
            end
        end
        function testLyapcholE(testCase) 
            for i=1:2
                if i<2
                    sys=loadSss('LF10'); % E is symmetric
                else
                    sys=loadSss('SpiralInductorPeec');
                end
                actR=lyapchol(sys);
                
                if sys.isDescriptor
                    actRes=norm(sys.A*(actR'*actR)*sys.E'+sys.E*(actR'*actR)*sys.A'+sys.B*sys.B');
                else
                    actRes=(sys.A*(actR'*actR)+(actR'*actR)*sys.A'+sys.B*sys.B');
                end
                
                actSolution={actRes};
                expSolution={zeros(size(actRes))};
                verification(testCase, actSolution, expSolution);
                
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-5,'AbsTol',1e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end