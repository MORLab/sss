classdef testSecondToFirst < sssTest
    
    methods (Test)  
        function testDifferentOptions(testCase)
            % test the differnt options for the F-matrix
            
            % option 'I'
            Opts.transf2nd = 'I';
            sys_I = sss('LF10',Opts);
            
            % option 'K'
            Opts.transf2nd = 'K';
            sys_K = sss('LF10',Opts);
            
            % option '-K'
            Opts.transf2nd = '-K';
            sys_minK = sss('LF10',Opts);
            
            % optoin 'alpha'
            Opts.transf2nd = 'alpha';
            sys_alpha = sss('LF10',Opts);
            
            % option scalar value
            Opts.transf2nd = 5;
            sys_scalar = sss('LF10',Opts);
            
            % option user defined matrix
            F = rand(18);
            F = (F+F');
            F = F + 20*eye(18);
            Opts.transf2nd = F;
            sys_matrix = sss('LF10',Opts);
            
            % create bode-plots
            [mag_I,phase_I,omega] = bode(sys_I);
            [mag_K,phase_K] = bode(sys_K,omega);
            [mag_minK,phase_minK] = bode(sys_minK,omega);
            [mag_alpha,phase_alpha] = bode(sys_alpha,omega);
            [mag_scalar,phase_scalar] = bode(sys_scalar,omega);
            [mag_matrix,phase_matrix] = bode(sys_matrix,omega);
            
            % check if the bode-plots are equal
            verifyEqual(testCase,mag_I,mag_K,'AbsTol',1e-8, ...
                'Verification failed for option "K" (magnitude)!');
            verifyEqual(testCase,phase_I,phase_K,'AbsTol',1e-8, ...
                'Verification failed for option "K" (phase)!');
            verifyEqual(testCase,mag_I,mag_minK,'AbsTol',1e-8, ...
                'Verification failed for option "-K" (magnitude)!');
            verifyEqual(testCase,phase_I,phase_minK,'AbsTol',1e-8, ...
                'Verification failed for option "-K" (phase)!');
            verifyEqual(testCase,mag_I,mag_alpha,'AbsTol',1e-8, ...
                'Verification failed for option "alpha" (magnitude)!');
            verifyEqual(testCase,phase_I,phase_alpha,'AbsTol',1e-8, ...
                'Verification failed for option "alpha" (phase)!');
            verifyEqual(testCase,mag_I,mag_scalar,'AbsTol',1e-8, ...
                'Verification failed for option "scalar" (magnitude)!');
            verifyEqual(testCase,phase_I,phase_scalar,'AbsTol',1e-8, ...
                'Verification failed for option "scalar" (phase)!');
            verifyEqual(testCase,mag_I,mag_matrix,'AbsTol',1e-8, ...
                'Verification failed for option "matrix" (magnitude)!');
            verifyEqual(testCase,phase_I,phase_matrix,'AbsTol',1e-8, ...
                'Verification failed for option "matrix" (phase)!');     
        end
        
        function testMatrixOutput(testCase)
            load('gyro')
            sys         =   second2first(M,1e-6*K,K,B,C);
            [A,B,C,D,E] =   second2first(M,1e-6*K,K,B,C);
            sys2 = sss(A,B,C,D,E);
            verifyClass(testCase, sys,'sss');
            verifyEqual(testCase, sys, sys2);        
        end
    end
end