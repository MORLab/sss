classdef testZpk < sssTest
    % testZpk - testing of zpk.m
   
    methods(Test)
        function testLM(testCase)
            
            systemsToTest{1} = loadSss('building');
            systemsToTest{2} = loadSss('heat-cont');
            systmesToTest{3} = loadSss('beam');
            
            for i=1:length(systemsToTest)
                sys=systemsToTest{i};
                sys_ss=ss(sys);

                if ~sys.isDae
                    k=20;
                    zpkData=zpk(sys,k,'lm');
                    zpkData_ss=zpk(sys_ss);
                    
                    %% verification of z
                    for j=1:sys.m
                        for l=1:sys.p      
                            if sys.isSiso
                                z=cell2mat(zpkData.z);
                                z_ss=cell2mat(zpkData_ss.z);
                            else
                                z=cell2mat(zpkData(j,l).z);
                                z_ss=cell2mat(zpkData_ss(j,l).z);
                            end
                            
                            z_ss = sort(z_ss,'descend');
                            z_ss = cplxpair(z_ss);
                            z = sort(z,'descend');
                            z = cplxpair(z);
                            
                            verification(testCase, z,z_ss(1:k));
                        end
                    end
                    
                    
                    %% verification of p
                    if ~sys.isSiso
                        p=cell2mat(zpkData(1,1).p);
                        p_ss=cell2mat(zpkData_ss(1,1).p);
                    else
                        p=cell2mat(zpkData.p);
                        p_ss=cell2mat(zpkData_ss.p);
                    end
                    p_ss = sort(p_ss,'descend');
                    p_ss = cplxpair(p_ss);
                    p = sort(p,'descend');
                    p = cplxpair(p);
                            
                    verification(testCase, p,p_ss(1:k));            
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-6,'AbsTol',1e-6,...
    'Difference between actual and expected exceeds relative tolerance');
end