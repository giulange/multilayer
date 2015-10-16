function compare_modules_versions(orig_mat, new_mat, type, name, range, pause_dur)
% compare_modules_versions(orig_mat, new_mat, type)
% 
% DESCRIPTION
%   This function compares simulations using different codes.
% 
% INPUT
%   orig_mat:       Mat file storing variables computed using the original
%                   model by Antonio Coppola (if type==1), or a model with
%                   oldest modifications in the code (if type==2).
%                   [mandatory]
% 
%   new_mat:        Mat file storing variables computed using the newest
%                   code implementation (both cases type==1 or type==2).
%                   [mandatory]
% 
%   type:           A flag to distinguish between the use of original
%                   Coppola's implementation (type==1) or another code (not
%                   the newest but newer than Coppola's one, with type==2).
%                   [mandatory]
% 
%   name:           If you set the name of the plot it will show you only
%                   that one! You can use one of the followings:
%                       °head
%                       °no3
%                       °nh4
%                   [optional]
% 
%   range:          Set the interval within which plot the simulation.
%                   [optional]
% 
%   pause_dur:      Duration in seconds of pause to use for any single plot
%                   of the simulation to be displayed. Default value is
%                   0.25 seconds.
%                   [optional]
% 
% EXAMPLE
%   compare_modules_versions('original/whole_sim_OLD_ver.mat', 'sim_pro-wt.mat', 1)
%   compare_modules_versions('original/whole_sim_OLD_ver.mat', 'sim_new-hsurf.mat', 1, 'no3')

%% def
figHW       = [ 350,600;    ...% head :: compare two versions
                0,0;        ...% ?
              ];
%% pre
% pause duration
if nargin < 6
    pause_dur = 0.25;
end
% Labels:
[~,Olabel,~] = fileparts(orig_mat);
[~,Nlabel,~] = fileparts(new_mat);
Olabel = ['REF{ ',Olabel,' }'];
Nlabel = ['NEW{ ',Nlabel,' }'];
%% PLOT HEAD PRESSURE :: COMPARE the two versions
% if ~exist('h1','var')
    % load from simulation using "original" version of multilayer:
%     load original/h1_orig.mat h1 time
% end
% if ~exist('O','var')
    % load from simulation using "new" version of multilayer:
%     load whole_sim.mat O P
%     load whole_sim_NEW_ver.mat O P
% end
if nargin < 4, name = 'head'; end
if type==1 && strcmp(name,'head')
    if type==1 && strcmp(name,'head')
        load(orig_mat,'h1','time')
    elseif type==2
        load(orig_mat,'O','P')
        time = P.time;
        h1 = O.h22;
    end
    load(new_mat,'O','P')

    Nnodes = size(O.h22,1);
    if size(h1,1)~=Nnodes, error('The number of nodes is not the same!'), end

    if nargin < 5 || isempty(range)
        dS = 0;
        dE = floor( max(P.time(1:P.j-1)) );
    else
        dS = max(0,range(1));
        dE = min(floor( max(P.time(1:P.j-1)) ),range(2));
    end

    Maxis   = [0,Nnodes,-180,2];

    for ii = dS:dE

        time_j = ii;

        % find original:
        [~,Fo] = min(abs(time-time_j));
        % find new:
        [~,Fn] = min(abs(P.time-time_j));

        figure(23),whitebg('k'),clf
        set(gcf,'Position',[ 1200, 100, figHW(1,2), figHW(1,1) ])

        if isfield(P,'sdate')
            timestamp = datestr(P.sdate+time_j,'yyyy-mmm-dd');
        else
            timestamp = num2str(time_j);
        end

        plot( 1:Nnodes,h1(:,Fo),'-or',1:Nnodes,O.h22(:,Fn),'-xg' )
        xlabel('Compartments','FontWeight','b','FontSize',14)
        ylabel('Head Pressure [cm]','FontWeight','b','FontSize',14)
        title(  {'\fontsize{18}Compare multilayer'; ...
                ['\fontsize{14}orig-mat vs new-mat version'];...
                ['\fontsize{15}\color{magenta}time = ',timestamp] }, ...
                'Interpreter','Tex')
        Caxis	= axis;

        Caxis(3) = min(Caxis(3),Maxis(3));
        Caxis(4) = max(Caxis(4),Maxis(4));
        axis( Caxis )

        hL = legend(Olabel,Nlabel);
        set(hL,'FontWeight','b','FontSize',8,'Interpreter','none', 'Location','SouthOutside')

        text(10,Caxis(3)+15, sprintf('time(%3d) = %9.4f',Fo,time(Fo)),'Color','red',...
             'BackGround',[.7,.7,.7])
        text(60,Caxis(3)+15, sprintf('time(%3d) = %9.4f',Fn,P.time(Fn)),'Color','green',...
             'BackGround',[.3,.3,.3])

        pause( pause_dur )
    end
end
%% PLOT HEAD PRESSURE :: new VERSION

% THIS CAN BE DONE INDIPENDENTLY FROM value in type!
% but I decide to skip:
if nargin < 4, name = 'head'; end
if type==3 && strcmp(name,'head')% I use type==3 to skip this

% if ~exist('O','var')
%     % load from simulation using "new" version of multilayer:
% %     load whole_sim.mat O P
%     load whole_sim_NEW_ver.mat O P
% end
load(new_mat,'O','P')

Nnodes = size(O.h22,1);

if nargin < 5
    dS = 0;
    dE = floor( max(P.time(1:P.j-1)) );
else
    dS = max(0,range(1));
    dE = min(floor( max(P.time(1:P.j-1)) ),range(2));
end

Maxis   = [0,Nnodes,-130,2];

for ii = dS:dE
    
    time_j = ii;

    % find new:
    [~,Fn] = min(abs(P.time-time_j));

    figure(24),whitebg('k'),clf
    plot( 1:Nnodes,O.h22(:,Fn),'-xg' )
    xlabel('Compartments','FontWeight','b','FontSize',14)
    ylabel('Head Pressure [cm]','FontWeight','b','FontSize',14)
    title(  {'\fontsize{18}Compare multilayer'; ...
            ['\fontsize{14}orig-mat vs new-mat version'];...
             ['\fontsize{15}\color{magenta}time = ',num2str(time_j+1)] },...
             'Interpreter','Tex')
    Caxis	= axis;
    
    Caxis(3) = min(Caxis(3),Maxis(3));
    Caxis(4) = max(Caxis(4),Maxis(4));
    axis( Caxis )
    
    hL = legend(Olabel,Nlabel);
    set(hL,'FontWeight','b','FontSize',8,'Interpreter','none', 'Location','SouthOutside')
    
    text(60,Caxis(3)+15, sprintf('time(%3d) = %9.4f',Fn,P.time(Fn)),'Color','green',...
         'BackGround',[.3,.3,.3])

    pause( pause_dur )
end
end
%% PLOT HEAD PRESSURE :: new VERSION :: before vs after modifications
if nargin < 4, name = 'head'; end
if type == 2 && strcmp(name,'head')
% load BEFORE:
% load whole_sim.mat O P
load(orig_mat,'O','P')
Ob = O;
Pb = P;
clear O P
% load AFTER:
% load whole_sim_NEW_ver.mat O P
% load whole_sim_temp.mat O P
load(new_mat,'O','P')
Oa = O;
Pa = P;
clear O P

Nnodes = size(Ob.h22,1);
if size(Oa.h22,1)~=Nnodes, error('The number of nodes is not the same!'), end

if nargin < 5
    dS = 0;
    dE = floor( max(Pb.time(1:Pb.j-1)) );
else
    dS = max(0,range(1));
    dE = min(floor( max(Pb.time(1:Pb.j-1)) ),range(2));
end

Maxis   = [0,Nnodes,-130,2];

for ii = dS:dE
    
    time_j = ii;

    % find BEFORE:
    [~,Fo] = min(abs(Pb.time-time_j));
    % find AFTER:
    [~,Fn] = min(abs(Pa.time-time_j));

    figure(25),whitebg('k'),clf
    plot( 1:Nnodes,Ob.h22(:,Fo),'-or',1:Nnodes,Oa.h22(:,Fn),'-xg' )
    xlabel('Compartments','FontWeight','b','FontSize',14)
    ylabel('Head Pressure [cm]','FontWeight','b','FontSize',14)
    title(  {'\fontsize{18}Compare multilayer'; ...
            ['\fontsize{14}orig-mat vs new-mat version'];...
             ['\fontsize{15}\color{magenta}time = ',num2str(time_j+1)] },...
             'Interpreter','Tex')
    Caxis	= axis;
    
    Caxis(3) = min(Caxis(3),Maxis(3));
    Caxis(4) = max(Caxis(4),Maxis(4));
    axis( Caxis )
    
    hL = legend(Olabel,Nlabel);
    set(hL,'FontWeight','b','FontSize',8,'Interpreter','none', 'Location','SouthOutside')
    
    text(10,Caxis(3)+15, sprintf('time(%3d) = %9.4f',Fo,Pb.time(Fo)),'Color','red',...
         'BackGround',[.7,.7,.7])
    text(60,Caxis(3)+15, sprintf('time(%3d) = %9.4f',Fn,Pa.time(Fn)),'Color','green',...
         'BackGround',[.3,.3,.3])

    pause( pause_dur )
end
end
%% PLOT NO3 :: COMPARE the two versions
% if ~exist('DNO','var')
%     % load from simulation using "original" version of multilayer:
%     load original/h1_orig.mat DNO time
% end
% if ~exist('O','var')
%     % load from simulation using "new" version of multilayer:
%     load whole_sim_NEW_ver.mat O P
% end

sl = 2;
if nargin < 4, name = 'no3'; end
if strcmp(name,'no3')
    if type==1
        load(orig_mat,'DNO','time')
        DNO = DNO(3:end,3:end);
    elseif type==2
        load(orig_mat,'O','P')
        DNO = O.C2(:,:,1,sl);
        time = P.time;
        Nnodes = size(O.h22,1);
        clear O P
    end

    load(new_mat,'O','P')
    if size(O.h22,1)~=Nnodes, error('The number of nodes is not the same!'), end

    if nargin < 5 || isempty(range)
        dS = 0;
        dE = floor( max(P.time(1:P.j-1)) );
    else
        dS = max(0,range(1));
        dE = min(floor( max(P.time(1:P.j-1)) ),range(2));
    end

    Maxis           = [0,Nnodes,0,1e-05];
    for ii = dS:dE

        time_j      = ii;

        % find original:
        if type==1
            Fo      = 1+ii;
        else
            [~,Fo]	= min(abs(time-time_j));
        end

        % find new:
        [~,Fn]      = min(abs(P.time-time_j));

        figure(26),whitebg('k'),clf
        plot( 1:Nnodes,DNO(:,Fo),'-or',1:Nnodes,O.C2(:,Fn,1,sl),'-xg' )
        xlabel('Compartments','FontWeight','b','FontSize',14)
        ylabel('NO_3^{-} [g / cm^2]','FontWeight','b','FontSize',14)
        title(  {'\fontsize{18}Compare multilayer';...
                ['\fontsize{14}orig-mat vs new-mat version'];...
                 ['\fontsize{15}\color{magenta}time = ',num2str(time_j)] },...
                 'Interpreter','Tex')

        Caxis       = axis;
        Caxis(3)    = min(Caxis(3),Maxis(3));
        Caxis(4)    = max(Caxis(4),Maxis(4));
        axis( Caxis )

        hL = legend(Olabel,Nlabel);
        set(hL,'FontWeight','b','FontSize',8,'Interpreter','none', 'Location','SouthOutside')

        text(80,Caxis(4)/3, sprintf('time(%3d) = %9.4f',Fo,ii),'Color','red',...
             'BackGround',[.7,.7,.7])
        text(80,Caxis(4)/10, sprintf('time(%3d) = %9.4f',Fn,P.time(Fn)),'Color','green',...
             'BackGround',[.3,.3,.3])

        pause( pause_dur )
    end
end
%% PLOT NH4 :: COMPARE the two versions
% if ~exist('DNH','var')
%     % load from simulation using "original" version of multilayer:
%     load original/h1_orig.mat DNH time
% end
% if ~exist('O','var')
%     % load from simulation using "new" version of multilayer:
%     load whole_sim_NEW_ver.mat O P
% end

sl=1;
if nargin < 4, name = 'nh4'; end
if strcmp(name,'nh4')
    if type==1
        load(orig_mat,'DNH','time')
        DNH = DNH(3:end,3:end);
    elseif type==2
        load(orig_mat,'O','P')
        DNH = O.C2(:,:,1,sl);
        time = P.time;
        Nnodes = size(O.h22,1);
        clear O P
    end

    load(new_mat,'O','P')
    if size(O.h22,1)~=Nnodes, error('The number of nodes is not the same!'), end

    if nargin < 5 || isempty(range)
        dS = 0;
        dE = floor( max(P.time(1:P.j-1)) );
    else
        dS = max(0,range(1));
        dE = min(floor( max(P.time(1:P.j-1)) ),range(2));
    end

    Maxis           = [0,Nnodes,0,1e-05];

    for ii = dS:dE

        time_j      = ii;

        % find original:
        if type==1
            Fo      = 1+ii;
        else
            [~,Fo]	= min(abs(time-time_j));
        end

        % find new:
        [~,Fn]      = min(abs(P.time-time_j));

        figure(26),whitebg('k'),clf
        plot( 1:Nnodes,DNH(:,Fo),'-or',1:Nnodes,O.C2(:,Fn,1,sl),'-xg' )
        xlabel('Compartments','FontWeight','b','FontSize',14)
        ylabel({'NH_4^{+} [g / cm^2]'},'FontWeight','b','FontSize',14)
        title(  {'\fontsize{18}Compare multilayer';...
                ['\fontsize{14}orig-mat vs new-mat version'];...
                 ['\fontsize{15}\color{magenta}time = ',num2str(time_j)] },...
                 'Interpreter','Tex')

        Caxis       = axis;
        Caxis(3)    = min(Caxis(3),Maxis(3));
        Caxis(4)    = max(Caxis(4),Maxis(4));
        axis( Caxis )

        hL = legend(Olabel,Nlabel);
        set(hL,'FontWeight','b','FontSize',8,'Interpreter','none', 'Location','SouthOutside')

        text(80,Caxis(4)/3, sprintf('time(%3d) = %9.4f',Fo,ii),'Color','red',...
             'BackGround',[.7,.7,.7])
        text(80,Caxis(4)/10, sprintf('time(%3d) = %9.4f',Fn,P.time(Fn)),'Color','green',...
             'BackGround',[.3,.3,.3])

        pause( pause_dur )
    end
end
%% end
return