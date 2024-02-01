%%%%
% Plotscript_Lambda - Generate figures of the lambda values
% throughout a BDCA run from the 'Logs' files resulting from the 
% 'Scaling Run' script. This is how those figures for the paper were
% generated. 
% 
% The parameters for the script: 
%
% Change the loaded file to select which 'Logs' file to generate figures
% for
% 
% select which cell(s) in the 'Logs' files to generate plots for
%%%%

% particular run to generate the plot from 
run = 1;

numruns = length([LogsDCA.time]);

%%
markersize = 8;
marker     = '.';

linemark = ':'; 

% Create the lambda plots
   
    %BDCA adaptive
    lambda = vertcat(LogsBDCA_adapt(run).logs.lambdaiter);
    taushifts = cumsum([LogsBDCA_adapt(run).logs.dcaiter]);
    
    figure;
    semilogy(lambda,marker,"MarkerSize",markersize);
    hold on 
    xline(taushifts(1:end-1),linemark);
    hold off
    axis padded
    xlabel('Total Iterations','Interpreter','latex');
    ylabel('$\lambda$','Interpreter','latex');
    title("BDCA Adaptive $\lambda$ values" ,'Interpreter','latex');
    fontsize(14,'points');