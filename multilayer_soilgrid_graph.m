function multilayer_soilgrid_graph(Pnodes, Wzint)
% multilayer_soilgrid_graph(Pnodes, Wzint)

%% pars
offset = [-10,+10,-20,+40];

%% main
scrsz = get(0,'ScreenSize');
fig_size = [scrsz(1)+20, scrsz(2), scrsz(3)*0.30, scrsz(4)*0.95];

Wzint = [0,Wzint];
figure(7)
clf
whitebg('k')
set(gcf, 'OuterPosition',fig_size)
hold on

% soil layers limits:
for ii = 1:length(Wzint)
    X = +[0,20];
    Y = -[Wzint(ii),Wzint(ii)];
    line(X,Y)
    line([X(2),X(2)+4],Y,'LineStyle','--','Color',[0.6,0.6,0.6])
    text(X(1),Y(1),num2str(Wzint(ii)),'HorizontalAlignment','right','VerticalAlignment','cap')
    if ii>1
        text(mean(X),-mean(Wzint(ii-1:ii)),['soil layer ',num2str(ii-1)],'HorizontalAlignment','center','VerticalAlignment','cap')
    end
end
% vertical lines of soil profile:
line([X(1),X(1)],-[Wzint(1),Wzint(end)])
line([X(2),X(2)],-[Wzint(1),Wzint(end)])
text(X(1),offset(4)/2,{'Depth';'[cm]'},'HorizontalAlignment','right','VerticalAlignment','bottom')

% soil grid:
for ii=1:length(Pnodes.z)-1
    text(X(2)+1,-Pnodes.z(ii),'.','Color','w','FontSize',6)
end
text(X(2)+1,-Pnodes.z(end),'.','Color','y','FontSize',6)

DIF = diff(Pnodes.soillayer(1:end-1));
DIF(find(DIF)-1) = 1;
DIF = [1;DIF];
DIF(end:end+1) = 1;
DIF = find(DIF);
for ii = 1:length(DIF)
    text(X(2)+1.5,-Pnodes.z(DIF(ii)),num2str(DIF(ii)),'Color','w','FontSize',6,'VerticalAlignment','cap')
    text(X(2)+3.5,-Pnodes.z(DIF(ii)),num2str(Pnodes.z(DIF(ii))),'Color','y','FontSize',6,'VerticalAlignment','cap')
end
text(X(2)+1.5,offset(4)/2,'nodes','HorizontalAlignment','left','VerticalAlignment','cap','rotation',90)
text(X(2)+3.5,offset(4)/2,{'Z';'[cm]'},'HorizontalAlignment','left','VerticalAlignment','bottom','Color','y')

axis([ X(1), X(2), -Wzint(end), Wzint(1) ]+offset)
hold off
axis off

%% end
return