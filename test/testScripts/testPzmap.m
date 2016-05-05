classdef testPzmap < sssTest
    % testPzmap - testing of pzmap.m
    %
    % Description:
    %   The function pzmap.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case the following inputs/outputs:
    %    pzmap(sys),[p,z]=pzmap(sys)
    %
    %
    % ------------------------------------------------------------------
    %   This file is part of sssMOR, a Sparse State Space, Model Order
    %   Reduction and System Analysis Toolbox developed at the Institute
    %   of Automatic Control, Technische Universitaet Muenchen.
    %   For updates and further information please visit www.rt.mw.tum.de
    %   For any suggestions, submission and/or bug reports, mail us at
    %                     -> sssMOR@rt.mw.tum.de <-
    % ------------------------------------------------------------------
    % Authors:      Alessandro Castagnotto, Maria Cruz Varona
    %               Jorge Luiz Moreira Silva
    % Last Change:  05 Nov 2015
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
   
    methods(Test)
        function test1(testCase)
            for i=1:length(testCase.sysCell)
                sysSparse=testCase.sysCell{i};
                sys=ss(sysSparse);
                [actP,actZ]=pzmap(sysSparse);
                
                % sort with abs values (cplxpair sorts with real values 
                % -> test of peec fails due to sorting)
                tbl=table(abs(actP),actP);
                tbl=sortrows(tbl);
                actP=tbl.Var1;
                tbl=table(abs(actZ),actZ);
                tbl=sortrows(tbl);
                actZ=tbl.Var1;
                actPZ={actP,actZ};
                
                % call without output arguments
                pzmap(sysSparse);
                
                % built-in
                [expP,expZ]=pzmap(sys);
                tbl=table(abs(expP),expP);
                tbl=sortrows(tbl);
                expP=tbl.Var1;
                tbl=table(abs(expZ),expZ);
                tbl=sortrows(tbl);
                expZ=tbl.Var1;
                expPZ={expP,expZ};
                verification(testCase, actPZ, expPZ);
                close all
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-6,'AbsTol',1e-6,...
    'Difference between actual and expected exceeds relative tolerance');
end