% Script for checking the effect of lambdaskip size on the timing results
clc;
clear; 

N = 10; % Run a test N times
skip_min = 0; 
skip_max = 500; 
skip_num = 100; 


[LogsDCA,SkipLogs,parameters] = skipping_check(N,skip_min,skip_max,skip_num);

filename = CreateUniqueFileName('profile/SkippingResults');
save(filename,'LogsDCA','SkipLogs','parameters');

skipping_plotting;