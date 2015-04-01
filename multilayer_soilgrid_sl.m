function Pnodes = multilayer_soilgrid_sl(Wsl,zint,Pnodes)
% Pnodes = multilayer_soilgrid_ud(Wsl,zint,Pnodes)
% 
% DESCRIPTION
%   This function creates a soil grid with node spacing defined by user
%   using the sub-layers method (see option "2" in W.sg.type).
%   The criterion used to create the grid is the following, assuming for
%   the first sub-layer (SoilLay=1) of the first layer (SubLay=1) 3 nodes
%   (nNodes=3) of 1 cm (hNode=1) with a total of 3 cm (hSubLay=3):
%   
%       0.0 -----   top-border (soil/atmosphere boundary)
%       0.5   *     1st node
%       1.0 -----   bottom-border (above node) / top-border (below node)
%       1.5   *     2nd node
%       2.0 -----   bottom-border (above node) / top-border (below node)
%       2.5   *     3rd node
%       3.0 -----   bottom-border (above node)
% 
% INPUT
%   Wsl:        It is the W.sg.sublayers variable passed as input to the
%               THIS function. It is made by the following 5 columns:
%                   -SoilLay:   The number of soil layer at which one or
%                               more sub-layers can be defined.
%                   -SubLay:    The global number of sublayer in soil
%                               profile.
%                   -hSubLay:   The thickness of current sublayer, which can
%                               be discretized into one or more nodes.
%                   -hNode:     The thickness of each node constituting the
%                               sublayer.
%                   -nNodes:    The number of regularly spaced nodes that
%                               discretise the sublayer.
% 
%   zint:       The depth of the lower border of soil layers (original
%               variable is W.zint).
% 
%   Pnodes:     The P.nodes variable is all passed in input as NaN. It is
%               used to store all the information to build the soil grid
%               used during the simulation in the multilayer program.
% 
% OUTPUT
%   Pnodes:     The P.nodes variable is returned filled with information
%               required to build the soil grid used during the multilayer
%               simulation.

%% pre-elaboration
iL                          = Wsl(:,1); % 'SoilLay'
iS                          = Wsl(:,2); % 'SubLay'
iHs                         = Wsl(:,3); % 'hSubLay'
iHn                         = Wsl(:,4); % 'hNode'
iNn                         = Wsl(:,5); % 'nNodes'
%% initialization
% Pnodes.num                  = NaN(  iL(end)+0,1 );
% Pnodes.thickness            = NaN(  iL(end)+0,1 );
% Pnodes.cumsum               = NaN(  iL(end)+1,1 );
% Pnodes.soillayer            = NaN( sum(iNn)+1,1 );
% Pnodes.z                    = NaN( sum(iNn)+1,1 );
% Pnodes.dz                   = NaN( sum(iNn)+1,1 );
%% check point
% The sum of sublayers thickness must equal the lower boundary depth of the
% lower soil horizon:
if zint(end) ~= sum( iHs )
    error('Error defininig hSubLay in the "sublayers" option!')
end
% The sequential number of soil sub-layers must continuous:
if sum( diff( iS )>1 )
    error('Error defininig SubLay in the "sublayers" option!')
end
% other checks?
% ...
%% main
for ii=1:iL(end)
    Pnodes.num(ii)          = sum(iNn(iL==ii));
    
    % It cannot be defined, because in the "sublayers" method nodes have
    % diffent thicknesses within the same soil layer:
    % Pnodes.thickness(ii)    = ?;
end

Pnodes.cumsum               = [0;cumsum(Pnodes.num)];

% ***TOP & INTERMEDIATE nodes within the soil grid***
% It counts the current number of node with soil profile
k                           = 0;
D                           = 0;
for ii=1:iS(end)
    for jj=1:iNn(ii)
        k                   = k + 1;
        D                   = D + iHn(ii);
        Pnodes.soillayer(k) = iL(ii);
        Pnodes.z(k)         = D - iHn(ii)*0.5;
        Pnodes.dz(k)        = iHn(ii);
    end
end
% ***BOTTOM node***
Pnodes.soillayer(end)       = NaN; % it already was NaN!
Pnodes.z(end)               = D + iHn(end)*0.5;
Pnodes.dz(end)              = iHn(end);
%% end
return