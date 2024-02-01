%%%%
% result - Test DCA, BDCA and BDCA adaptive for a constrained
% clustering problem on the continental data set. This generates 3 
% variables 'LogsBDCA', 'LogsBDCA_adapt', and 'LogsDCA' which contain
% detailed information on all of the solves. These Log files are then saved
% to a file, the profile directory, which can then be be plotted by
% 'plotscript'
% 
%
% See Example 7.4 in the paper.
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
