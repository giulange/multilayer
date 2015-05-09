% clean & clear
clear,clc

% reset profiling tool:
profile clear

% activate profiling:
profile on -history

try
    
    % run program to profile:
    multilayer_run.m

    % stop profiling
    profile off

    % record profiling info
    p = profile('info');

    % explore profiled program:
    profile viewer
    
    % save profile data:
    profsave(p,'profile_results')

catch exception

    % stop profiling
    profile off

end