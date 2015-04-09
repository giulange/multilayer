%% Initialization
run multilayer_init.m
%% MONTECARLO --START--
for mm=1:M.nvp
    progress_bar(mm,M.nvp,'multilayer prog')
%     hwb = waitbar(0,['Please wait... [MCS=',num2str(mm),']']);
%% Preallocation
    run multilayer_preallocation
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
%% SOIL GRID GEOMETRY
    switch W.sg.type
        case 1 % regular soil grid
            P.nodes             = multilayer_soilgrid(P.nz,W.zint,P.nodes);
        case 2 % sl --> sub-layers
            P.nodes             = multilayer_soilgrid_sl(W.sg.sublayers,W.zint,P.nodes);
    end
    % plot: (should be optional)
%     multilayer_soilgrid_graph(P.nodes,W.zint);
%% Create printing matrices -- check with Antonio
    % useless?
    W.dt                = W.dtin;
    
    % *** -- put in initialization? -- ***
    % Inizializzazione contatori e tempo di simulazione
    % (mettiamo i contatori in struct array specifica come "C"!)
    P.j                 = 1;
    P.time(1)           = W.dt;
    P.kk                = 1;
    P.pp                = 1;
    P.L                 = 0;
    P.LL                = 0;
    P.CC                = 0;
    P.T                 = 0;
%     P.TT                = 0;
    P.rnf               = 0;
    P.ktec              = 1;  
%% Calcolo ET potenziale con E + T -- è già scritto dentro while su P.kk
% ETp0, Ep and Tp and Droot will be defined again later at P.kk=1:
%     if and(W.itopvar==1,W.iveg==1)
%         P.ETp0          = V.Kc(1)*V.ETr(1);
%         % legge di Beer
%         P.Ep            = P.ETp0*exp(-V.extf*V.LAI(1));
%         P.Tp            = P.ETp0-P.Ep;
%         P.Droot         = V.Droot(1);
%     end
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
        P.sh.ifc(inode)        = W.ifc(    P.nodes.soillayer(inode) );
        % Calcolo del teta alla field capacity da utilizzare nel
        % modulo soluti per il calcolo del fattore di riduzione del
        % coefficiente di mineralizzazione.
        P.sh.tetafc(inode)     = fnteta( W.hfc, P.sh, inode );
    end
%% TIME simulation loop
    while P.time(P.j)<W.tmax
%         waitbar( P.j / P.Nj )
%% controllo sui tempi di stampa -- check:: position within while? correct?
        if abs(P.time(P.j)-W.tp(P.pp))<0.0001
            P.TT                = 1;
        elseif P.time(P.j)>W.tp(P.pp)
            P.time(P.j)         = W.tp(P.pp);
            W.dt                = W.tp(P.pp)-P.time(P.j-1);
            P.TT                = 1;
        else
            P.TT                = 0;
        end
%% adattamento T sim al T stampa -- check:: position within while? correct?
        if P.time(P.j) >= W.tp(P.pp)
            P.TT                = 1;
            P.time(P.j)         = W.tp(P.pp);
            W.dt                = W.tp(P.pp);
        elseif W.tp(P.pp)-P.time(P.j) < W.dtmin && W.tp(P.pp)-P.time(P.j)>=0
            P.TT                = 1;
            P.time(P.j)         = W.tp(P.pp);
            W.dt                = W.tp(P.pp);
        end
%% controllo sui flussi al contorno superiore -- ATTENTION: dovrebbe stare in water transport (solve tq!)
if W.itopvar==1 && P.L==0 % L ==> tutto ok nel j precedente
    P.tq                    = B.top.thqstar(P.kk+1);

    if W.itbc==0                    % flux
        W.qsurf             = B.top.hqstar(P.kk);
    elseif and(W.itbc==1,P.rnf==0)  % potential
        W.hsurf             = B.top.hqstar(P.kk);
    end

    if W.ibotvar==1
        if W.ibbc==0                % flux
            W.qbot          = B.bot.hqstar(P.kk);
        elseif W.ibbc==1            % potentials
            W.hbot          = B.bot.hqstar(P.kk);
        end
    end

    if and(W.isol==2,W.iCtopvar==1)
        P.Cinput(1)         = B.Ctop.Cstar.NH.FR(P.kk);
        P.Cinput(2)         = B.Ctop.Cstar.NO.FR(P.kk);
    end

    if W.iveg==1
        P.ETp0              = V.Kc(P.kk)*V.ETr(P.kk);
        P.Ep                = P.ETp0*exp(-V.extf*V.LAI(P.kk));
        P.Tp                = P.ETp0-P.Ep;
        P.Droot             = V.Droot(P.kk);
    end
end        
%% controllo tempi top boundary -- check:: position within while? correct?
        if W.itopvar==1
            if abs(P.time(P.j)-P.tq)<0.00001
                P.T             = 1;
                P.L             = 0;
            elseif P.time(P.j)>P.tq
                P.time(P.j)     = P.tq;
                W.dt            = P.tq-P.time(P.j-1);
                P.T             = 1;
                P.L             = 0;
            else
                P.T             = 0;
                P.L             = 1;
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
        run multilayer_transport_water_us.m
%% ?? define cell ??
        P.jstar                     = P.j;
        if P.LL==0
%% Solutes Transport
            if W.isol==2
            % update intial concentrations
                if P.j==1% first iteration
                    % [W.dz,1:2,P.Nj]
                    P.C1(:,1)       = S.CDE.Cin.NH;
                    P.C1(:,2)       = S.CDE.Cin.NO;
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
%% restore qsurf & hsurf -- check with Antonio!!
% Serve a ripristinare il valore di W.qsurf che potrebbe essere stato
% cambiato nella routine per l'evaporazione W.qsurf=P.Ep oppure
% W.qsurf=P.Emax.
% Lo stesso vale per W.hsurf che potrebbe essere stato cambiato nella
% routine per W.itbc=1 W.hsurf=(3*P.h1(1)-P.h1(2))/2 oppure
% W.hsurf=W.hsurfmax.
% Se questi valori non venissero ripristinati, nel giro successivo del
% while i controlli sui flussi in superficie (per W.itbc=0) o sui
% potenziali in superficie (per W.itbc=1) verrebbro effettuati non
% utilizzando i valori letti nel file di input ma su quelli calcolati nella
% routine.
            if and(W.itopvar==1,W.itbc==0)
                W.qsurf=B.top.hqstar(P.kk);
            elseif and(and(W.itopvar==1,W.itbc==1),P.rnf==0)
                W.hsurf=B.top.hqstar(P.kk);
            elseif and(W.itbc==1,P.rnf==1)
                W.qsurf=B.top.hqstar(P.kk);
            end
%             W.hqsurf     = B.top.hqstar(P.kk);
%             W.hsurf     = B.top.hqstar(P.kk);
%% runoff [& runon]
            O.runoff(1,P.j,mm)   = W.qsurf-O.fluxsurf(1,P.j,mm);
%% print potentials [useless ??]
            % aggiorna il contatore per il tempo di stampa (ci serve
            % ancora??)
            if P.TT==1
                P.pp        = P.pp+1;                
            end
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
%   Ovviamente questa verifica va fatta solo se si � entrati in W.itbc=1
%   partendo da W.itbc=0.
            if P.T==1               % flag counter top-bound
                P.kk=P.kk+1;        % contatore top-bound & Ctop-bund  
                 if W.itopvar==1    % che vuol dire?
                     % se il flusso al nuovo kk > kk-1 ho ancora runoff
                     if and(P.rnf==1,abs(B.top.hqstar(P.kk))<abs(B.top.hqstar(P.kk-1)))
                        P.rnf=0;    % potrebbe non avere piu' senso 
                        W.itbc=0;
                     end
                 end
            end
%% update time of simulation
            P.j         = P.j+1;
            P.time(P.j) = P.time(P.j-1)+W.dt;
        end% if P.LL=0    
%% PLOT -- tmp
        figure(88),whitebg('k')
        hold on,subplot(411),plot([P.sh.tetafc,P.teta])
        legend('tetafc','teta'),title(sprintf('j = %4d',P.j-1)), hold off;
        hold on,subplot(412),plot([P.sink]), legend('sink'),hold off;
        hold on,subplot(413),plot([P.h1,P.h2]), legend('h1','h2'),hold off;
        hold on,subplot(414),plot([P.cap,P.kond]), legend('cap','kond'),hold off;
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