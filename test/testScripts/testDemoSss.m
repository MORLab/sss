classdef testDemoSss < sssTest
    
    methods (Test)  
        function mainFunctionality(testCase)
            Opts.pause=false;
            Opts.test=true;
            sss_gettingStarted(Opts); % in sss
        end
    end
end