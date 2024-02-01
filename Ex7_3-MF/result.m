%%%%
% result - Test DCA, BDCA and BDCA adaptive for a constrained
% clustering problem on the EIL76 data set. This generates 3 
% variables 'LogsBDCA', 'LogsBDCA_adapt', and 'LogsDCA' which contain
% detailed information on all of the solves. These Log files are then saved
% to a file, the profile directory, which can then be be plotted by
% 'plotscript'
%
% See Example 7.3 in the paper.
%
% The parameters for the script
%
% N  - Number of times to rerun for a different random initial point.
%
% Output: 
%
% LogsBDCA, LogsBDCA_adapt, LogsDCA   - 
% 
%       Log files containing detailed information about all of the
%       runs. A plotting script 'plotscript' can use these log files to 
%       generate figures of the same kind used in the paper. 
%%%%

clc;
clear; 

N = 100; % Run a test N times

[LogsDCA,LogsBDCA_adapt,parameters] = loopmain(N);

filename = CreateUniqueFileName('profile/Results');
save(filename,'LogsDCA','LogsBDCA_adapt','parameters');

%plotting 
plotscript;


function[FileName]=CreateUniqueFileName(FileName)
[fPath, fName, fExt] = fileparts(FileName);
if isempty(fExt)  % No '.mat' in FileName
  fExt     = '.mat';
  FileName = fullfile(fPath, [fName, fExt]);
end
if exist(FileName,'file')
    [fPath, fName, fExt] = fileparts(FileName);
    fDir = dir(fullfile(fPath, [fName,' (*)', fExt]));
    if isempty(fDir)
        FileName = fullfile(fPath, [fName,' (1)', fExt]);
    else
        pattern = "(" + digitsPattern + ")" + fExt;
        hasThePattern = endsWith(extractfield(fDir,'name'),pattern);
        Extracted = extract(extractfield(fDir(hasThePattern),'name'),pattern);
        num = max(cell2mat(cellfun(@(C) textscan(C,'(%d)') , Extracted,'UniformOutput',true)));
        num = num+1;
        FileName = fullfile(fPath, [fName,' (',num2str(num),')', fExt]);
    end
end
end
