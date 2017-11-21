classdef testLyapchol < sssTest
    methods(Test)
        function testLyapchol1(testCase) 
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    [S,R]=lyapchol(sys);
                    verification(testCase, sys, S, R);
                end
            end
        end
        function testLyapcholOpts(testCase)
            Opts.method='adi';
            Opts.q=120;
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && sys.n>300
                    [S,R]=lyapchol(sys,Opts);
                    verification(testCase, sys, S, R);
                end
            end
        end
        function testLyapcholE(testCase) 
            for i=1:3
                if i<2
                    sys=sss('LF10'); % E is symmetric
                    Opts.method = 'auto';
                else
                    sys=sss('SpiralInductorPeec');
                end
                if i == 2
                    Opts.method='adi';
                elseif i == 3
                    Opts.method='crksm';
                    Opts.lowrank = 1;
                end    
                [S,R]=lyapchol(sys,Opts);               
                verification(testCase, sys, S, R);      
            end
        end
        function testMatrices(testCase) 
            sys  = sss('SpiralInductorPeec');
            
            % first test adi-method
            Opts.method = 'adi';
            [S,R]=lyapchol(sys,Opts);  
            verification(testCase, sys, S, R);  

            % compare to matrix call (should be the same)
            S2 = sssFunc.lyapchol(sys.A,sys.B,sys.E);
            R2 = sssFunc.lyapchol(sys.A',sys.C',sys.E');
            verification(testCase, sys, S2, R2); 
            verifyEqual(testCase,S,S2,'AbsTol',1e-3);
            verifyEqual(testCase,R,R2,'AbsTol',1e-3);

            % only 2 matrices (E=I)
            sysI = sys; sysI.E = speye(sysI.n);
            S = sssFunc.lyapchol(sysI.A, sysI.B);
            R = sssFunc.lyapchol(sysI.A',sysI.C');
            verification(testCase, sysI, S, R); 
            clear S R S2 R2
            
            % secomd test crksm-method
            Opts.method = 'crksm';
            Opts.lowrank = 1;
            [S,R]=lyapchol(sys,Opts);  
             verification(testCase, sys, S, R);  
             
              % compare to matrix call (should be the same)
            S2 = sssFunc.lyapchol(sys.A,sys.B,sys.E,Opts);
            R2 = sssFunc.lyapchol(sys.A',sys.C',sys.E',Opts);
            verification(testCase, sys, S2, R2); 
            verifyEqual(testCase,S,S2,'AbsTol',1e-3);
            verifyEqual(testCase,R,R2,'AbsTol',1e-3);

            % only 2 matrices (E=I)
            sysI = sys; sysI.E = speye(sysI.n);
            S = sssFunc.lyapchol(sysI.A, sysI.B,Opts);
            R = sssFunc.lyapchol(sysI.A',sysI.C',Opts);
            verification(testCase, sysI, S, R); 
           
        end
    end
end

function [] = verification(testCase,sys,  S, R)
    X = S*S'; Y = R*R';
    tol = 1e-3;
    verifyLessThan(testCase, ...
        norm(sys.A*(X)*sys.E' + sys.E*(X)*sys.A' + sys.B*sys.B'),...
        tol);
    verifyLessThan(testCase, ...
        norm(sys.A'*(Y)*sys.E + sys.E'*(Y)*sys.A + sys.C'*sys.C),...
        tol);

end