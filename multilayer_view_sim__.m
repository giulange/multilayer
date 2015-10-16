function multilayer_view_sim__( sim_matFile )
% multilayer_view_sim__( sim_matFile )
% 
% DESCRIPTION
%   This function plots the pressure head, the nitrate concentration and
%   the ammonium concentration of a full/partial simulation using the
%   multilayer program.
% 
% INPUT
%   sim_matFile:    The saved mat file in which the required simulation by
%                   multilayer is stored.

%% config
backColor   = [0.4,0.4,0.4];
frontColor  = {'yellow', 'red', 'green'};
FIL_psh     = '/home/giuliano/work/Projects/LIFE_Project/swap/soil_hydr.mat';
%% load
[~,Fname,~] = fileparts(sim_matFile);
load( sim_matFile, 'P', 'O' )
%% pre
timeday     = NaN(P.tidx,1);
figHW       = [850,450];
%% plot LOCATION
if isfield(P,'X')
    % Soil Hydraulic Functions
    load( FIL_psh )
    flERROR     = true;

    % --- find current polygon from COORDINATES
    SOIL_UNIT   = [];
    for ii = 1:numel(SHAPE)
        [IN,ON] = inpolygon( P.X, P.Y, SHAPE(ii).X, SHAPE(ii).Y );
        if IN || ON
            SOIL_UNIT = ii;
            flERROR = false;
            break
        end
    end
    if ~strcmpi( SOIL(ii).SIGLA_UC, P.SIGLA_UC ), flERROR = true; end

    % --- find current polygon from soil unit class LABEL
    FNS = find( strcmpi(P.SIGLA_UC, {SOIL.SIGLA_UC}) );
    if ~isempty(SOIL_UNIT)
        flERROR = false;
    end 

    if flERROR
        error('coordinates :: check!')
    end

    figure(2),clf,hold on
    %                   [left,bottom,width,height]
    subplot('Position', [0.0, 0.1,   0.9,  0.9])
            mapshow(SHAPE,'DisplayType','polygon','FaceColor','w')
            mapshow(SHAPE(SOIL_UNIT),'DisplayType','polygon','FaceColor','g')
            for ii = 1:length(FNS)
                if FNS(ii)==SOIL_UNIT, continue, end
                mapshow(SHAPE(FNS(ii)),'DisplayType','polygon','FaceColor','y')
            end
            mapshow(P.X,P.Y,'DisplayType','point')
        axis equal, axis off
    %                    [left,bottom,width,height]
        axes('Position', [0.7, 0.0,   0.3,  0.3])
        title( ['\fontsize{8}\color{green}SOIL UNIT :: ', SOIL(SOIL_UNIT).SIGLA_UC,' (N=',num2str(numel(FNS)),')'] )
            mapshow(SHAPE(SOIL_UNIT),'DisplayType','polygon','FaceColor','g')
            mapshow(P.X,P.Y,'DisplayType','point', 'Marker','+','MarkerEdgeColor', 'r','MarkerFaceColor','r','MarkerSize',8, 'LineWidth',2)
        axis equal, axis off
    hold off
end
%% main
for iji = 1:P.tidx-1
    if isempty( find(P.time==iji,1) )
%         timeday(iji)=[];
        continue
    end
    timeday(iji) = find(P.time==iji,1);
end

% here I have NaNs where the crisp date is not present but there is the
% crisp date +/- the tolerance. I should implement the search for the
% non-crisp date, instead of removing/skipping NaNs as I do now!
% timeday(isnan(timeday))=[];

figure(17),clf,whitebg('k')
% [left, bottom, width, height]
set(gcf,'Position',[ 20, 20, figHW(2), figHW(1) ])
for iji=1:length(timeday)
    % Skip NaNs:
    if isnan(timeday(iji)), continue, end
    figure(17)
    subplot(311),hold on,
        for i_iji = 1:iji-1
            % Skip NaNs:
            if isnan(timeday(i_iji)), continue, end
            plot(O.h22(:,timeday(i_iji)),'color',backColor)
        end
        plot(O.h22(:,timeday(iji)),'color',frontColor{1}),hold off
        ylabel('Head Pressure [cm]','FontWeight','b','FontSize',14)

    if isfield(P,'sdate')
        timestamp = datestr(P.sdate+iji,'yyyy-mmm-dd');
    else
        timestamp = num2str(iji);
    end
    if isfield(P,'ElapsedTime_of_Simulation')
        ElapsedTime = ['\fontsize{14}\color{white}Elapsed Time: ',sprintf('%2d[min] %2d[sec]',floor(P.ElapsedTime_of_Simulation/60),round(mod(P.ElapsedTime_of_Simulation,60)))];
    else
        ElapsedTime = '';
    end
    title( {
             '\fontsize{18}View Simulation'; ...
            ['\fontsize{14}\color{green}File :: ',strrep(Fname, '_','-')];...
            ElapsedTime; ...
            ['\fontsize{15}\color{magenta}time = ',timestamp]
           }, 'Interpreter','Tex')

    subplot(313),hold on,
        for i_iji = 1:iji-1
            % Skip NaNs:
            if isnan(timeday(i_iji)), continue, end
            plot(O.C2(:,timeday(i_iji),1,2),'color',backColor)
        end
        plot(O.C2(:,timeday(iji),1,2),'color',frontColor{2}),hold off
%         xlabel('Nodes','FontWeight','b','FontSize',14)
        ylabel('NO^{-}_{3} [g/cm^3]','FontWeight','b','FontSize',14)
    subplot(312),hold on,
        for i_iji = 1:iji-1
            % Skip NaNs:
            if isnan(timeday(i_iji)), continue, end
            plot(O.C2(:,timeday(i_iji),1,1),'color',backColor)
        end
        plot(O.C2(:,timeday(iji),1,1),'color',frontColor{3}),hold off
        xlabel('Nodes','FontWeight','b','FontSize',14)
        ylabel('NH^{+}_{4} [g/cm^3]','FontWeight','b','FontSize',14)
    pause(0.5)
end
%% exit
return