% run multilayer_assign_vars % calling --> multilayer_conf.m
% NOW I use directly multilayer_conf.m for simplicity:

%% pars
fpf = @(TXT) fprintf('%s\n', TXT);

%% prog

fpf('Loading prog vars...')
run multilayer_conf_admin.m

fpf('Loading user vars...')
run multilayer_conf.m

fpf('Checking vars...')
run multilayer_checkpars.m

% *****TEMPORARY*****
fpf('Input data...  --->(temporary condition waiting for an exhaustive input framework!!)')
run multilayer_data_tmp.m
% *****TEMPORARY*****

fpf('Running prog...')
run multilayer_prog.m

% fpf('Displaying results...')
% run multilayer_graph.m

%% clean
clear fpf