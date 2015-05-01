%% Initialization
run multilayer_init.m
%% MONTECARLO --START--
for mm=1:M.nvp
    progress_bar(mm,M.nvp,'multilayer prog')
%     hwb = waitbar(0,['Please wait... [MCS=',num2str(mm),']']);
%% Preallocation
    run multilayer_preallocation.m
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
%% Create printing matrices -- check with Antonio
    % useful:
    P.dt                = W.dtin;
    P.time(1)           = P.dt;
    P.rnf               = 0;
    P.ktec              = 1;
    P.tidx_jm1          = floor(P.time(1))+1;
    figure(98),clf
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
        % time element to extract from time-dependent parameter-vectors:
        P.tidx              = floor(P.time(P.j))+1;
        P.L                 = P.tidx > P.tidx_jm1;
        % set the check upon the print times
        if P.L
            if P.time(P.j)-P.tidx >= W.tptole && P.tidx-P.time(P.j-1) >= W.tptole
                % Set the proper dt increment for current simulation to
                % print at the time defined by user plus the tolerance.
                dtprevious  = P.dt;
                P.time(P.j) = P.tidx; %P.time(P.j) -dtprevious +P.dt;
                P.dt        = P.tidx-P.time(P.j-1);
            end
        end
%         waitbar( P.j / P.Nj )
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
            run multilayer_transport_water_us.m
        else        
            % execute a different water transport model(s):
            %  -for instance only for saturated:
            run multilayer_transport_water_s.m
            %  -or only unsaturated:
            run multilayer_transport_water_u.m
        end
%% ?? define cell ??
        P.jstar                 = P.j;
%% Solutes Transport
        if W.isol==2
        % update intial concentrations
            if P.j==1% first iteration
                % [W.dz,1:2,P.Nj]
                P.C1(:,1)       = P.CDECinNH;
                P.C1(:,2)       = P.CDECinNO;
            elseif P.j>P.jstar% se avanza...
                % OLD % P.C1(:,:,P.j) = P.C2(:,:,P.j-1);
                P.C1            = squeeze(O.C2(:,P.j-1,mm,:));   % UGUALI??
            elseif P.j==P.jstar% se non puo' avanzare...
                % OLD % P.C1(:,:,P.j) = P.C1star(:,:,P.j);
                P.C1            = squeeze(O.C2(:,P.j-1,mm,:));   % UGUALI??
            end
            %[O,P] = multilayer_transport_solute_N_ade( P, W, S, B, O, mm );
            [O.C2(:,P.j,mm,:),P]= multilayer_transport_solute_N_ade(P,W,S.CDE,...
                                  B.Ctop,P.C1,O.fluxsurf(1,P.j,mm),...
                                  O.fluxbot(1,P.j,mm) );
        end
%% runoff [& runon]
        O.runoff(1,P.j,mm)      = W.qsurf-O.fluxsurf(1,P.j,mm);
%% controllo di W.itbc
%   Questo controllo va fatto solo se W.itopvar=1.
%   Se W.itopvar=0, una volta avvenuto il cambio da W.itbc=0 a W.itbc=1,
%   non e' piu' necessario tornare ad W.itbc=0 perche' il flusso non cambia
%   fino a fine simulazione.
%   Se il nuovo B.top.hqstar(P.kk) e' maggiore del vecchio
%   (B.top.hqstar(P.kk-1) non e' necessario cambiare W.itbc da 1 a 0.
%   Se invece e' minore allora bisogna di nuovo verificare se il nuovo
%   B.top.hqstar sia maggiore di P.fluxsurf_max, ed allora occorre cambiare
%   W.itbc da 1 a 0.
%   Ovviamente questa verifica va fatta solo se si è entrati in W.itbc=1
%   partendo da W.itbc=0.
        if P.L==1               % flag counter top-bound
            % se il flusso al nuovo P.tidx > P.tidx-1 ho ancora runoff
            if and(P.rnf==1,abs(B.top.hqstar(P.tidx))<abs(B.top.hqstar(P.tidx_jm1)))
                P.rnf=0;        % potrebbe non avere piu' senso 
                W.itbc=0;
            end
        end
%% update time of simulation
        P.j         = P.j+1;
        P.time(P.j) = P.time(P.j-1)+P.dt;
        P.tidx_jm1  = P.tidx;
%% PLOT -- tmp
%         figure(88),whitebg('k')
%         hold on,subplot(411),plot([P.sh.tetafc,P.teta])
%         legend('tetafc','teta'),title(sprintf('j = %4d',P.j-1)), hold off;
%         hold on,subplot(412),plot([P.sink]), legend('sink'),hold off;
%         hold on,subplot(413),plot([P.h1(1:10),O.h22(1:10,P.j)]), legend('h1','h2'),hold off;
%         hold on,subplot(414),plot([P.cap,P.kond]), legend('cap','kond'),hold off;
        figure(98),whitebg('k'),
        subplot(211),hold on,plot(O.h22(:,P.j-1)),hold off,legend('O.h22 cumulative');
        title(sprintf('time(%4d) = %10.3f',P.j-1,P.time(P.j-1)),'FontSize',14,'FontWeight','b')
        subplot(212),plot(O.h22(:,P.j-1)),legend('O.h22');
    end% P.time(P.j)<W.tmax
%% MONTECARLO --END--
end% mm=1:M.nvp
% close(hwb) % close waitbar
%% SAVE -- incomplete [set which times and nodes to print!!]
if W.MTCL == 0 || M.nvp == 1
    multilayer_save( O, proj, P.jstar )
elseif W.MTCL == 1
    multilayer_save_mcs( O, proj )
end