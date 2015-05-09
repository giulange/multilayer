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
    P.rnf               = 0;
    P.ktec              = 1;
    P.tidx_jm1          = floor(P.time(1))+1;
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
        % 'ifc' is the same along the whole soil profile, hence it must be
        % transformed into a scalar!!
        P.sh.ifc(inode)        = W.ifc(    P.nodes.soillayer(inode) );
        % Calcolo del teta alla field capacity da utilizzare nel
        % modulo soluti per il calcolo del fattore di riduzione del
        % coefficiente di mineralizzazione.
        P.sh.tetafc(inode)     = fnteta( W.hfc, P.sh, inode );
    end
%% TIME simulation loop
    while P.time(P.j)<W.tmax
%% update tidx & dt
        % time element to extract from time-dependent parameter-vectors:
        P.tidx              = floor(P.time(P.j))+1;
%       use mod() to simulate the last piece of day till P.tidx, without
%       updating input vars to the new DAY
        
        P.L                 = P.tidx > P.tidx_jm1;
        % set the check upon the print times
        if P.L
%             waitbar( P.j / P.Nj )

            % DELETE the following if to allow the printing time exactly at
            % P.tidx time!! Be aware that the printing at the exact time is
            % not ensured when the tridiagonal system cannot be solved at
            % current dt (and the program use a different dt
            if P.time(P.j)-P.tidx >= W.tptole && P.tidx-P.time(P.j-1) >= W.tptole
                % Set the proper dt increment for current simulation to
                % print at the time defined by user plus the tolerance.
                dtprevious  = P.dt;
                P.time(P.j) = P.tidx; %P.time(P.j) -dtprevious +P.dt;
                P.dt        = P.tidx-P.time(P.j-1);
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
            multilayer_transport_water_s
            %  -or only unsaturated:
            multilayer_transport_water_u
        end
%% Solutes Transport
        if W.isol==2
        % update intial concentrations
            if P.j==1% first iteration
                % [W.dz,1:2,P.Nj]
                P.C1(:,1)       = P.CDECinNH;
                P.C1(:,2)       = P.CDECinNO;
            else
                % OLD % P.C1(:,:,P.j) = P.C2(:,:,P.j-1);
                P.C1            = squeeze(O.C2(:,P.j-1,mm,:));
            end
            %[O,P] = multilayer_transport_solute_N_ade( P, W, S, B, O, mm );
            [O.C2(:,P.j,mm,:),P]= multilayer_transport_solute_N_ade(P,W,S.CDE,...
                                    B.Ctop,P.C1,O.fluxsurf(1,P.j,mm),...
                                    O.fluxbot(1,P.j,mm), O.h22(:,P.j,mm) );
        end
%% Update :: j, time, tidx
        P.j         = P.j+1;
        P.time(P.j) = P.time(P.j-1)+P.dt;
        P.tidx_jm1  = P.tidx;
    end% P.time(P.j)<W.tmax
%% MONTECARLO --END--
end% mm=1:M.nvp
% close(hwb) % close waitbar
%% SAVE -- incomplete [set which times and nodes to print!!]
if W.MTCL == 0 || M.nvp == 1
    multilayer_save( O, proj, P.j-1 )
elseif W.MTCL == 1
    multilayer_save_mcs( O, proj )
end