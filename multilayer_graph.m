%% plot results in time using the soil grid graph

%% Select solute:
% sl = 1; % ammonia NH4
sl = 2; % nitrate NO3

%% plot soil grid
multilayer_soilgrid_graph(P.nodes,W.zint);
%% transport processes
hold on
% other hydrological components of system:
color_in = 'white';
color_out= 'black';
% -pond(triangle):
o = [2.0,0];
fill( [0,1.0,-1.0,0]+o(1), [0,5,5,0]+o(1), 'g' ) % starter:bottom vertex
% -evaporation (arrow)
o = [5,0];
line( [0,0,1.2,0,-1.2,0]+o(1), [0,10,10,14,10,10]+o(2), 'Color',color_out ) % starter:bottom vertex
text(0+o(1),15+o(2),'Evapor','HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','n','FontSize',8,'Color',color_out)
text(0+o(1),24+o(2),num2str(0.001),'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','d','FontSize',8,'EdgeColor',color_out)
% -traspiration
o = [9,0];
line( [0,0,1.2,0,-1.2,0]+o(1), [0,10,10,14,10,10]+o(2), 'Color',color_out ) % starter:bottom vertex
text(0+o(1),15+o(2),'Transp','HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','n','FontSize',8,'Color',color_out)
text(0+o(1),24+o(2),num2str(0.002),'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','d','FontSize',8,'EdgeColor',color_out)
% -rain
o = [13,0];
line( [0,0,1.2,0,-1.2,0]+o(1), [14,4,4,0,4,4]+o(2), 'Color',color_in ) % starter:top vertex
text(0+o(1),15+o(2),'Rain','HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','n','FontSize',8,'Color',color_in)
text(0+o(1),24+o(2),num2str(0.003),'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','d','FontSize',8,'EdgeColor',color_in)
% -irrigation
o = [17,0];
line( [0,0,1.2,0,-1.2,0]+o(1), [14,4,4,0,4,4]+o(2), 'Color',color_in ) % starter:top vertex
text(0+o(1),15+o(2),'Irrigat','HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','n','FontSize',8,'Color',color_in)
text(0+o(1),24+o(2),num2str(0.004),'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','d','FontSize',8,'EdgeColor',color_in)
% -interception
% o = [10,25];
% line( [0,0,1.2,0,-1.2,0]+o(1), [14,4,4,0,4,4]+o(2), 'Color','w' ) % starter:top vertex

% -runoff
o = [28,0];
% line( [0,-4,-4,-6,-4,-4]+o(1), [0,0,3,0,-3,0]+o(2), 'Color','w' ) % starter:right vertex
line( [0,4,4,6,4,4]+o(1), [0,0,3,0,-3,0]+o(2), 'Color',color_out ) % starter:left vertex
text(+8.5+o(1),-4+o(2),'runoff','HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','n','FontSize',8,'Color',color_out)
text(+2+o(1),6+o(2),num2str(0.005),'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','d','FontSize',8,'EdgeColor',color_out)

% -runon
o = [-10,-0];
line( [0,4,4,6,4,4]+o(1), [0,0,3,0,-3,0]+o(2), 'Color',color_in ) % starter:left vertex
text(-2+o(1),-4+o(2),'runon','HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','n','FontSize',8,'Color',color_in)
text(+2+o(1),6+o(2),num2str(0.006),'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','d','FontSize',8,'EdgeColor',color_in)

% -infiltration/exfiltration

% -up/downward seepage

% -lateral flow

hold off
%% plot solute concentration simulation
if sl==1
    sol_lab = 'NH_h^{+}';
elseif sl==2
    sol_lab = 'NO_3^{-}';
end

[X,Y] = meshgrid([1,19],-(P.nodes.z(1:end-1)-P.nodes.dz(1:end-1)*0.5));
Csolute = squeeze( O.C2(:,1:P.jstar,1,sl) );

% Csolute = log(Csolute +0.0001);

mima = minmax( Csolute(:)' );
% cmap = colormap(jet(128));
% caxis( mima )
hold on
colorbar
for tj = 1:P.jstar
    title( {'{\fontsize{16}\color{yellow}NO_3^{-}} \fontsize{13}\color{green}', ...
        sprintf('Timestep: %4d',tj)} )
    Z = repmat( Csolute(:,tj), 1, 2);
    pcolor(X,Y,Z)
    pause(0.5)
end
hold off