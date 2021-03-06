% run multilayer_assign_vars % calling --> multilayer_conf.m
% NOW I use directly multilayer_conf.m for simplicity:

%% pars
fpf = @(TXT) fprintf('%s\n', TXT);

%% prog

% fpf('Static initialization...')
% run multilayer_preallocation_static.m

fpf('Loading prog vars...')
multilayer_conf_admin

fpf('Loading user vars...')
multilayer_conf

fpf('Checking vars...')
multilayer_check_and_load

% *****TEMPORARY*****
% % fpf('Input data...  --->(temporary condition waiting for an exhaustive input framework!!)')
% % multilayer_data_tmp
% *****TEMPORARY*****

fpf('Running prog...')
multilayer_prog

% SAVE -- incomplete [set which times and nodes to print!!]
fpf('Saving sim outputs...')
if W.save
    if W.MTCL == 0 || M.nvp == 1
        multilayer_save( O, proj, P.j-1, W.isol )
    elseif W.MTCL == 1
        multilayer_save_mcs( O, proj )
    end
end

if proj.video
    fpf('Displaying results...')
    multilayer_graph
end
%% clean
clear fpf