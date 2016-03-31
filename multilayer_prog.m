tic;
%% Pre-allocation
multilayer_init
%% MONTECARLO --START--
for mm=1:M.nvp
%     progress_bar(mm,M.nvp,'multilayer prog')
%     hwb = waitbar(0,['Please wait... [MCS=',num2str(mm),']']);
%% Preallocation
    multilayer_preallocation
%% Define conditions for current montecarlo setting
    if W.MTCL==1
        % ripristinare uso di M.nnc
        
        % Assign the stochastic values to the variables to be simulated
        % using Montecarlo:
        for nvars=1:length(M.list)
            % set current montecarlo condition:
            eval( [M.list{ nvars } ' = M.data(M.combinations(',num2str(mm),',:),',num2str(nvars),')'';'] )
        end
    end
%% Initialize :: dt, time, tidx, rnf, ktec 
    % useful:
    P.dt                = W.dtin;
    P.time(1)           = P.dt;
    P.tidx_jm1          = floor(P.time(1))+1;
    P.tidx              = floor(P.time(1))+1;
    P.L                 = P.tidx > P.tidx_jm1;
%     hFig = figure(98);clf
%     set(hFig,'Name','Potentials & time increments','Position',[800, 20, 500, 800]),whitebg('k'),
%% Hydraulic mapping :: soil layers --> soil grid nodes
    % Retention & conductivity at soil grid nodes from info at soil layers:
    for inode = 1:P.nz
        P.sh.dap(inode)        = W.dap(    P.nodes.soillayer(inode) );
        P.sh.tetas(inode)      = W.tetas(  P.nodes.soillayer(inode) );
        P.sh.tetar(inode)      = W.tetar(  P.nodes.soillayer(inode) );
        P.sh.alfrs(inode)      = W.alfrs(  P.nodes.soillayer(inode) );
        P.sh.fi(inode)         = W.fi(     P.nodes.soillayer(inode) );
        P.sh.alfvg(inode)      = W.alfvg(  P.nodes.soillayer(inode) );
        P.sh.en(inode)         = W.en(     P.nodes.soillayer(inode) );
        P.sh.alfvg2(inode)     = W.alfvg2( P.nodes.soillayer(inode) );
        P.sh.en2(inode)        = W.en2(    P.nodes.soillayer(inode) );
        P.sh.ifr(inode)        = W.ifr(    P.nodes.soillayer(inode) );
        P.sh.k0(inode)         = W.k0(     P.nodes.soillayer(inode) );
        P.sh.k0macr(inode)     = W.k0macr( P.nodes.soillayer(inode) );
        P.sh.bita(inode)       = W.bita(   P.nodes.soillayer(inode) );
        P.sh.bita2(inode)      = W.bita2(  P.nodes.soillayer(inode) );
        P.sh.h_enpr(inode)     = W.h_enpr( P.nodes.soillayer(inode) );
        % 'ifc' is the same along the whole soil profile, hence it must be
        % transformed into a scalar!!
        P.sh.ifc(inode)        = W.ifc(    P.nodes.soillayer(inode) );
        % Calcolo del teta alla field capacity da utilizzare nel
        % modulo soluti per il calcolo del fattore di riduzione del
        % coefficiente di mineralizzazione.
        %P.sh.tetafc(inode)     = fnteta( W.hfc, P.sh, inode );
        P.sh.tetafc(inode)     = multilayer_fnteta_vgm( W.hfc, P.sh, inode );
    end
%% Init SoilWater State rates/variables
    if W.wt_mod==1% New solver based on Newton-Raphson algorithm
        multilayer_soilwater_init
        if W.isol ~= 0
            multilayer_solute_init
        end
    end
%% TIME simulation loop
    fprintf( '%s\n', repmat('_',1,130) );
    fprintf('%14s, %7s %4s %11s %7s %11s %3s %5s %5s %6s %6s %8s %12s %12s\n', ...
        'date','P.time','P.j','dtwater','1/P.dt','dtsolute','p','rain','irri','roots','requ','stress','camm_day','cnit_day')
    while P.time(P.j)<W.tmax
%% CROP :: rotation & irrigation
        if W.iveg && P.flStartOfDay
% ROTATION   Schedule :: Load/Unload Crop(s)
            multilayer_rotation_schedule
% IRRIGATION Schedule :: flirri=[1,2,3]
            if V.flirri>0
                multilayer_irrigation
            end
        end
%% calcolo valori potenziale osmotico per il tempo di simulazione (ERROR)
        if W.iosm==1 && V.ifs>3
            P.IEC               = 0;        % boolean: store if ktec incremented
            %# ERROR: Antonio checks!! --> P.time(P.j)>EC.t(P.ktec)
            if P.time(P.j)>EC.t(P.ktec+1)
                P.ktec          = P.ktec+1;
                P.IEC           = 1;
            end
            P.ECstar            = ( P.time(P.j)-EC.t(P.ktec) )             /...
                                  ( EC.t(P.ktec+1)-EC.t(P.ktec) )          *...
                                  ( EC.data(:,P.ktec+1)-EC.data(:,P.ktec) )+...
                                    EC.data(:,P.ktec) ;
        end
%% Water Transport
        if true
            % execute saturated/unsaturated model:
            multilayer_transport_water_us
        else
            % execute a different water transport model(s):
            %  -for instance only for saturated:
            multilayer_transport_water_s %#ok<UNRCH>
            %  -or only unsaturated:
            multilayer_transport_water_u
        end
%% Solutes Transport
        if W.isol==2
% === OLD module ======================================
%         % update intial concentrations
%             if P.j==1% first iteration
%                 % [W.dz,1:2,P.Nj]
%                 P.C1(:,1)       = P.CDECinNH;
%                 P.C1(:,2)       = P.CDECinNO;
%             else
%                 % OLD % P.C1(:,:,P.j) = P.C2(:,:,P.j-1);
%                 P.C1            = squeeze(O.C2(:,P.j-1,mm,:));
%             end
%            %[O,P] = multilayer_transport_solute_N_ade( P, W, S, B, O, mm );
%            [O.C2(:,P.j,mm,:),P]= multilayer_transport_solute_N_ade(P,W,S.CDE,...
%                                    B.Ctop,P.C1,O.fluxsurf(1,P.j,mm),...
%                                    O.fluxbot(1,P.j,mm), O.h22(:,P.j,mm) );
%             solbal = zeros(1,2);
%             P.Ndtsolute(P.tidx) = P.dt;
% === NEW module ======================================
            multilayer_transport_solute
            O.C2(:,P.j,mm,:) = P.cml;% [g cm-3 water]
            P.solbal(:,P.j) = solbal;
        end
%% CROP :: growth model (inactive)
        % This module is put before timecontrol because P.flEndOfDay is
        % switched true before the last j of the current day is peformed!
        if P.flEndOfDay
            fprintf( '%14s, %7d %4d %11.6f %7d %11.6f %3d %5.2f %5.2f %6.2f %6.2f %8.2f %12.8f %12.8f\n', ...
                datestr(datenum(W.sdate) +P.tidx -1,'yyyy-mmm-dd'), ...
                round(P.time(P.j)),P.j,P.dt,round(1/P.dt),1/P.Ndtsolute(P.tidx),P.iter(1,P.j), ...
                B.top.rain(P.tidx), P.irri(P.tidx),...
                P.Droot(P.tidx), P.dstor(P.tidx), P.hrz_cm(P.tidx),solbal(1),solbal(2));
            if W.iveg, multilayer_cropgrowth, end
        end
%% Update :: { j, dt, time, tidx, flEndOfDay }
        multilayer_timecontrol
    end% P.time(P.j)<W.tmax
    fprintf('\n%14s, %7s %4s %11s %7s %3s %11s %5s %5s %6s %6s %8s %12s %12s\n', ...
        'date','P.time','P.j','dtwater','1/P.dt','p','dtsolute','rain','irri','roots','requ','stress','camm_day','cnit_day');
    fprintf( '%14s: %7d %4d %11.6f %7d %3d %11.6f %5.2f %5.2f %6.2f %6.2f %8.2f %12.8f %12.8f\n', ...
        'TOT', ...
        NaN,NaN,NaN,NaN,NaN,NaN, ...
        sum(B.top.rain), sum(P.irri),...
        P.Droot(P.tidx), nansum(P.dstor), max(P.hrz_cm), ...
        nansum(P.solbal(1,:)),nansum(P.solbal(2,:)) );
    fprintf( '%s\n', repmat('_',1,130) );
%% MONTECARLO --END--
end% mm=1:M.nvp
% close(hwb) % close waitbar
P.ElapsedTime_of_Simulation = toc;
%% prepare "O" for saving purpose
% time, tprint, tptolle
O.time      = P.time;
O.tprint    = P.tprint;
O.tptolle   = W.tptole;
O.z         = P.nodes.z(1:end-1);% last index is for bottom boundary!
% store teta
for ii = 1:size(O.h22,2)
    O.teta(:,ii) = multilayer_fnteta_vgm( O.h22(:,ii), P.sh, 1:P.nz );
end
