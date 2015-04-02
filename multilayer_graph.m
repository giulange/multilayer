%% plot results in time using the soil grid graph

%% Select solute:
% sl = 1; % ammonia NH4
sl = 2; % nitrate NO3

%% plot soil grid
multilayer_soilgrid_graph(P.nodes,W.zint);
% pause(2)
%% plot solute concentration simulation
[X,Y] = meshgrid([1,19],-[P.nodes.z(1:end-1)-P.nodes.dz(1:end-1)*0.5]);
Csolute = squeeze( O.C2(:,1:P.jstar,1,sl) );

% Csolute = log(Csolute +0.0001);

mima = minmax( Csolute(:)' );
% cmap = colormap(jet(128));
% caxis( mima )
hold on
% colorbar
for tj = 1:P.jstar
    title( {['{\fontsize{16}\color{yellow}NO_3^{-}} \fontsize{13}\color{green}'], ...
        sprintf('Timestep: %4d',tj)} )
    Z = repmat( log(Csolute(:,tj)), 1, 2);
    pcolor(X,Y,Z)
    pause(0.5)
end