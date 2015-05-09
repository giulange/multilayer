function runoff = multilayer_runoff(pond,pondmax,rsro,rsroexp,dt)
%% lateral drainage
% I assume no drainage:
swdra = 0;
% I should parameterize this item when lateral drainage will be implemented
% I set wls equal to zero (because lateral drainage is missing):
wls = 0;
% storage in sw-reservoir, per unit area
swst = 0;
% function swstlev(.,.): to be implemented! --> see functions.for
%% main

% --- if (POND is less than or equal to PONDMAX):
runoff = 0;

% --- condition ALWAYS satisfied!
if (pond-pondmax>0) && (swdra~=2)
    if rsro < 1.0d-03
        runoff = pond-pondmax;
    else
        % Eq. 4.2, page 71, SWAP-32:
        runoff = dt/rsro * (pond-pondmax)^rsroexp;
    end
    
% --- condition NEVER satisfied!
elseif swdra==2
    if pond>pondmax && pond>wls
        runoff = dt/rsro * (pond-max(pondmax,wls))^rsroexp;
    elseif pond<wls
        inun_max = swst - swstlev(pond,sttab);
        runoff = min(inun_max,wls-max(pond,pondmax));
    end
end
%% exit
return