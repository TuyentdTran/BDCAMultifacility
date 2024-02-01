%%%%
% Plotscript_Lambda - Generate figures of the lambda values
% from the NoSkipVsSkipBDCA files
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

numruns = length([LogsBDCA_adapt.time]);

%%
markersize = 8;
marker     = '.';

linemark = ':'; 

% Create the lambda plots
   
    %BDCA adaptive skipping
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
    title("BDCA Adaptive $\lambda$ Values With Skipping " ,'Interpreter','latex');
    fontsize(14,'points');


      %BDCA adaptive no skip
    lambda = vertcat(LogsBDCA_adapt_noskip(run).logs.lambdaiter);
    taushifts = cumsum([LogsBDCA_adapt_noskip(run).logs.dcaiter]);
    
    figure;
    semilogy(lambda,marker,"MarkerSize",markersize);
    hold on 
    xline(taushifts(1:end-1),linemark);
    hold off
    axis padded
    xlabel('Total Iterations','Interpreter','latex');
    ylabel('$\lambda$','Interpreter','latex');
    title("BDCA Adaptive $\lambda$ Values Without Skipping" ,'Interpreter','latex');
    fontsize(14,'points');