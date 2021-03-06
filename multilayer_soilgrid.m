function Pnodes = multilayer_soilgrid(numnodes,zint)
% Pnodes = multilayer_soilgrid( numnodes, zint )
% 

%% pre-elaboration
nlay                        = length(zint);
zdepth                      = zint - [0,zint(1:end-1)];
%% initialization
Pnodes.num                  = NaN( nlay+0,    1 );
Pnodes.thickness            = NaN( nlay+0,    1 );
Pnodes.cumsum               = NaN( nlay+1,    1 );
Pnodes.soillayer            = NaN( numnodes+1,1 );
Pnodes.z                    = NaN( numnodes+1,1 );
Pnodes.dz                   = NaN( numnodes+1,1 );
%% main
% Mono- and Multi- stratum:
for ii=1:nlay
    Pnodes.num(ii)         = round(numnodes * zdepth(ii) / zint(nlay));
    if ii==nlay
        if numnodes - sum(Pnodes.num)==0
            % Correct!!
        elseif numnodes - sum(Pnodes.num)==1
            % Adjust!!
            Pnodes.num(ii)  = Pnodes.num(ii)+1;
        else
            % Exit!!
            error('Wrong number of nodes!!')
        end
    end
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
Pnodes.soillayer(end+1)    = NaN; % use circshift instead!
Pnodes.z(end)              = Pnodes.z(end-1) + 2*(zint(nlay)-Pnodes.z(end-1));
Pnodes.dz(end)             = Pnodes.z(ii) - Pnodes.z(ii-1);
%% end
return