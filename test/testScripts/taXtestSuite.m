function taXtestSuite()
% This test suite is executing the test models located in the testcases
% directory of taX. It reports tests that fail.
clear;clc;close all;

% Go into tax root folder in order to be able to execute taXinit and git
cd(fileparts(mfilename('fullpath')));
taXinit

%% Execute regression tests
% Default maximum frequency
fMax = 101;
% Default tolerance
tol = 10^-10;
% Default frequencies to be tested
omega = [-10i+10 50i-3 101i+7 ]*2*pi;

% Switch to reference branch
!git checkout TestReference
!git submodule update
% Retrieve and Filter simulink models
path = mfilename('fullpath');
path = fileparts(path);
testCases = getTestcases(fullfile(path,filesep,'Testcases'));
testCases = testCases(~cellfun('isempty',strfind(testCases,'.slx')),:);

runTest(testCases,true,fMax,tol,omega)

% Switch to current master branch
!git checkout master
!git submodule update
runTest(testCases,false,fMax,tol,omega)

end

function runTest(testCases,regenerateReference,fMax,tol,omega)
% Iterate over all available testcases
for i = 1: size(testCases,1);
    modelPath = char(testCases(i));
    try
        sys = tax(modelPath,fMax);
    catch
        sys = sss();
    end
    
    [path, modelName, ~] = fileparts(modelPath); % get file name
    referenceName = [path,filesep,'reference_',modelName,'.mat'];
    
    if regenerateReference
        % Save result
        save(referenceName,'sys') %#ok<*UNRCH>
        disp('Regenerated Results');
    else
        passed = true;
        % Load reference
        reference = load(referenceName); %#ok<*UNRCH>
        
        %% Compare to reference
        % System Matrix
        if (det(sys.A)-det(reference.sys.A))^2 < tol
            disp('Passed det(A).')
        else
            warning('Failed to pass det(A).')
            disp(['Model:',modelPath])
            disp(['Reference:',referenceName])
            passed = false;
        end
        % Frequency response/transfer function
        if sum(sum(sum((freqresp(sys,omega)-freqresp(reference.sys,omega)).^2))) < tol
            disp('Passed freqresp.')
        else
            warning('Failed to pass freqresp.')
            disp(['Model:',modelPath])
            disp(['Reference:',referenceName])
            passed = false;
        end
        if passed
            % Remove reference
            delete(referenceName);
        end
    end
end
end

function fileList = getTestcases(dirName)
dirData = dir(dirName);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
fileList = fileList(~cellfun('isempty',strfind(fileList,'.slx')),:);
if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
        fileList,'UniformOutput',false);
end
end
