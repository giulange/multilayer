function Pnodes = multilayer_soilgrid(numnodes,zint, Pnodes)
% Pnodes = multilayer_soilgrid(numnodes,zint, P.nodes)

% SUBSTITUTIONS:
% W.nz      --> W.numnodes          % in future release!
% P.dztop   --> P.nodes.dz(1)       % done!
% P.dzbot   --> P.nodes.dz(W.nz+1)  % done!
% P.d_z     --> P.nodes.thickness   % done!
% P.istar   --> P.nodes.cumsum      % done!
% P.z       --> P.nodes.z           % done!
% P.dz      --> P.nodes.dz          % done!
% P.zmax    --> W.zint(W.nlay)      % not found (?!)

%% main
nlay                        = length(zint);
zdepth                      = zint - [0,zint(1:end-1)];

% Mono- and Multi- stratum:
for ii=1:nlay
    Pnodes.num(ii)         = round(numnodes * zdepth(ii) / zint(nlay));
    Pnodes.thickness(ii)   = zdepth(ii) / Pnodes.num(ii);
end

Pnodes.cumsum              = [0;cumsum(Pnodes.num)];

% TOP node of the soil grid:
Pnodes.z(1)                = Pnodes.thickness(1)/2;
Pnodes.dz(1)               = Pnodes.thickness(1)/2;
% INTERMEDIATE nodes of the soil grid:
for jj=1:nlay
   for ii=Pnodes.cumsum(jj)+2:Pnodes.cumsum(jj+1)+1
       Pnodes.soillayer(ii)= jj;
       Pnodes.z(ii)        = Pnodes.z(ii-1) + Pnodes.thickness(jj);
       Pnodes.dz(ii)       = Pnodes.z(ii) - Pnodes.z(ii-1);
   end
end
% BOTTOM node of the soil grid:
Pnodes.soillayer(1)        = [];
Pnodes.z(end)              = Pnodes.z(end-1) + 2*(zint(nlay)-Pnodes.z(end-1));
Pnodes.dz(end)             = Pnodes.z(ii) - Pnodes.z(ii-1);
%% end
return