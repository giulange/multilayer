% multilayer_transport_water_su
% 
% DESCRIPTION
%   This script simulates the transport of water in unsaturated/saturated
%   soils at a time step specified by the main "multilayer" routine.
%   It solves the Richard's equation at each compartment of the soil
%   system.
%   ...top boundary implementation
%   ...bottom boundary implementation
%   ...more hydrological details on what is performed...
%   It can handle a fixed and a variable compartment height.
%   ...other characteristics of the model implemented.

%% init
% Newton-Raphson/Thomas was broken by an exit condition:
nr_breaked  = false;
% Newton-Raphson/Thomas convergence not possible at this timestep:
fl_noconv   = false;
% Newton-Raphson/Thomas number of non-convergences:
n_noconv    = 0;
%% Module for water transport
switch W.wt_mod
    case 0
        % run only if P.timestep raise of 1 integer step!
%         if P.L==1
%             run multilayer_boundtop_th.m
%         end
        while ~nr_breaked && ~fl_noconv
            run multilayer_boundtop_th.m
            run multilayer_wtus_thomas.m
        end
    case 1
        while ~nr_breaked && ~fl_noconv
            run multilayer_wtus_newtonraphson.m
        end
    
end
%% clean
clear alllllllll