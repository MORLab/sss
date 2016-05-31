classdef testZpk < sssTest
    % testZpk - testing of zpk.m
   
    methods(Test)
        function testLM(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};

                if ~sys.isDae
                    k=20;
                    [p,zAll]=zpk(sys,k);
                    if ~sys.isSiso
                        p=p{1,1};
                    end
                    
                    %% verification of z
                    for k=1:sys.m
                        for l=1:sys.p      
                            if sys.isSiso
                                z=zAll;
                            else
                                z=zAll{k,l};
                            end

                            % compare results of eig and eigs
                            zEig=eig(full([sys.A,sys.B(:,l);sys.C(k,:),sys.D(k,l)]),[full(sys.E),zeros(sys.n,1);zeros(1,sys.n),0]);

                            % remove zeros at infinity
                            zEig=zEig(abs(real(zEig))<1e6);

                            tbl=table(-abs(z),z);
                            tbl=sortrows(tbl);
                            z=tbl.z;
                            tbl=table(-abs(zEig),zEig);
                            tbl=sortrows(tbl);
                            zEig=tbl.zEig;
                            
                            %disp([z,zEig(1:size(z,1)),abs(z),abs(zEig(1:size(z,1)))]);

                            % remove elements not contained in z from zEig (but all
                            % elements of z must be in zEig)
                            j=1;
                            nRemoved=0;
                            
                            while(j<size(z,1))
                                if  abs(zEig(j,1)-z(j,1))>abs(z(j,1))*1e-4
                                    zEig(j)=[];
                                    j=1;
                                    nRemoved=nRemoved+1;
                                else
                                    j=j+1;
                                end
                                if j==size(zEig,1)
                                    break
                                end
                            end

                            disp(['Elements not contained in z: ',num2str(nRemoved)]);

                            %remove single complex element (for comparison)
                            if abs(imag(sum(z)))>1e-12
                               z(abs(imag(z)-imag(sum(z)))<1e-16)=[];
                            end

                            verifyLessThanOrEqual(testCase, nRemoved,10);
                            verification(testCase, z,zEig(1:size(z,1)));
                        end
                    end
                    
                    
                    %% verification of p
                    
                    % compare results of eig and eigs
                    pEig=eig(full(sys.A),full(sys.E));

                    % remove poles at infinity
                    pEig=pEig(abs(real(pEig))<1e6);

                    tbl=table(-abs(p),p);
                    tbl=sortrows(tbl);
                    p=tbl.p;
                    tbl=table(-abs(pEig),pEig);
                    tbl=sortrows(tbl);
                    pEig=tbl.pEig;

                    %remove single complex element (for comparison)
                    if abs(imag(sum(p)))>1e-12
                       p(abs(imag(p)-imag(sum(p)))<1e-16)=[];
                    end
                    
                    % make sure sorting is correct
                    p=cplxpair(p);
                    pEig=cplxpair(pEig(1:size(p,1)));
                    
                    verification(testCase, p,pEig(1:size(p,1)));
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-6,'AbsTol',1e-6,...
    'Difference between actual and expected exceeds relative tolerance');
end