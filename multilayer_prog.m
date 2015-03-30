%% initialization

% We need a standard:
%   possible indexes are --> { node, time, montecarlo, solute? }
%   These are the rules:
%   ---to-be-saved---
%   (1) node, (2) time, (3) montecarlo, (4) solute;
%   (1) node, (2) time, (3) montecarlo;
%   (1) 1,    (2) time, (3) montecarlo;
%   ---temporary---
%   (1) node, (2) time;
%   (1) 1,    (2) time;
%   (1) node, (2) 1;
%   other...

% ****TO BE SAVED****
%   > O.C2          --> solutes concentrations
%                       [nodes x times x montecarlo x solutes]
%   NOTE: I have to substitute C2 with a variable specific for any solute
%         we would like to implement in the model.
O.C2                            = NaN( W.nz, P.Nj, M.nvp, 2 );
%   > O.h22         --> flussi ai nodi intermedi
%                       [nodes x times x montecarlo]
O.h22                           = NaN( W.nz, P.Nj, M.nvp );
%   > O.fluxsurf    --> flusso al contorno superiore
%                       [1 x times x montecarlo]
O.fluxsurf                      = NaN( 1, P.Nj, M.nvp );
%   > O.fluxbot     --> flusso al contorno inferiore
%                       [1 x times x montecarlo]
O.fluxbot                       = NaN( 1, P.Nj, M.nvp );
%   > O.runoff      --> runoff
%                       [1 x times x montecarlo]
O.runoff                        = NaN( 1, P.Nj, M.nvp );

% ****TEMPORARY****
% --nodes
P.nodes.num                     = NaN(W.nlay+0,1);
P.nodes.thickness               = NaN(W.nlay+0,1);
P.nodes.cumsum                  = NaN(W.nlay+1,1);
P.nodes.soillayer               = NaN(W.nz+0,1);
P.nodes.z                       = NaN(W.nz+1,1);
P.nodes.dz                      = NaN(W.nz+1,1);
% --soil-grid
P.sh.dap                        = NaN(W.nz,1);
P.sh.tetas                      = NaN(W.nz,1);
P.sh.tetar                      = NaN(W.nz,1);
P.sh.alfrs                      = NaN(W.nz,1);
P.sh.fi                         = NaN(W.nz,1);
P.sh.alfvg                      = NaN(W.nz,1);
P.sh.en                         = NaN(W.nz,1);
P.sh.alfvg2                     = NaN(W.nz,1);
P.sh.en2                        = NaN(W.nz,1);
P.sh.ifr                        = NaN(W.nz,1);
P.sh.k0                         = NaN(W.nz,1);
P.sh.k0macr                     = NaN(W.nz,1);
P.sh.bita                       = NaN(W.nz,1);
P.sh.bita2                      = NaN(W.nz,1);
P.sh.ifc                        = NaN(W.nz,1);
P.sh.tetafc                     = NaN(W.nz,1);
% --scalars:
% ---counters:
P.j                             = NaN;
P.jstar                         = NaN;
P.k                             = NaN;
P.SS                            = NaN;
P.niter                         = NaN;
% ---others:
P.dpt                           = NaN;
P.op                            = NaN;
P.teta_hsurf                    = NaN;
P.Emax                          = NaN;
P.teta_hbot                     = NaN;
% --vectors:
P.h2                            = NaN(W.nz,1);
% ---time:
P.time                          = NaN(1,P.Nj);
P.km_max                        = NaN(1,P.Nj);
P.fluxsurf_max                  = NaN(1,P.Nj);
P.km                            = NaN(1,P.Nj);
% ---others:
P.Cinput                        = NaN(2,1);
% --arrays: % check with Antonio --> I would DELETE the P.Nj dimension!!
P.C1                            = NaN(W.nz,P.Nj,2);
P.h1                            = NaN(W.nz,2,P.Nj); % **check with Antonio
P.h1star                        = NaN(W.nz,2,P.Nj); % **check with Antonio
P.ECstar                        = NaN(W.nz,P.Nj);
P.teta                          = NaN(W.nz,P.Nj);
P.kond                          = NaN(W.nz,P.Nj);
P.cap                           = NaN(W.nz,P.Nj);
P.sink                          = NaN(W.nz,P.Nj);
P.kp                            = NaN(W.nz,P.Nj);
P.flux                          = NaN(W.nz,P.Nj);
%% MONTECARLO --START--
for mm=1:M.nvp
    progress_bar(mm,M.nvp,'multilayer prog')
%% define conditions (e.g. soil layers) in current montecarlo setting
    if W.MTCL==1
        % ripristinare uso di M.nnc
        
        % Assign the stochastic values to the variables to be simulated
        % using Montecarlo:
        for nvars=1:length(M.list)
            % ANTONIO: IS THAT FINE? IS WHAT WE NEED?
            eval( [M.list{ nvars } ' = M.data(M.combinations(',num2str(mm),',:),',num2str(nvars),')'';'] )
        end
    end
%% Definition of soil grid geometry
    P.nodes = multilayer_soilgrid(W.nz,W.zint,P.nodes);
    
    % You should call it as:
    %     P = multilayer_soilgrid(W.nz,W.zint,P);
    % when you'll set W.zint without the 0 !!
%% Create printing matrices -- put in initialization --
    
    W.dt                = W.dtin;
    
    % inizializzazione contatori e tempo di simulazione
    P.j                 = 1; % mettiamo i contatori in struct array specifica come "C"
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
%% Calcolo ET potenziale con E + T
    if and(W.itopvar==1,W.iveg==1)
        P.ETp0          = V.Kc(1)*V.ETr(1);
        % legge di Beer
        P.Ep            = P.ETp0*exp(-V.extf*V.LAI(1));
        P.Tp            = P.ETp0-P.Ep;
        P.Droot         = V.Droot(1);
    end
%% retention & conductivity at soil grid nodes
    for inode = 1:W.nz
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
%% inizio simulazione
    while P.time(P.j)<W.tmax
%% adattamento T sim al T stampa
        if P.time(P.j) >= W.tp(P.pp)
            P.TT = 1;
            P.time(P.j)     = W.tp(P.pp);
            W.dt            = W.tp(P.pp);
        elseif W.tp(P.pp)-P.time(P.j) < W.dtmin && W.tp(P.pp)-P.time(P.j)>=0
            P.TT = 1;
            P.time(P.j)     = W.tp(P.pp);
            W.dt            = W.tp(P.pp);
        end    
%% lettura h iniziali o aggiornamento h
        P.jstar                 = P.j;
        
        % lettura h iniziali e aggiornamento h
        if P.j==1
            P.h1(:,P.j)         = W.hin;
        elseif P.j>P.jstar
            P.h1(:,P.j)         = P.h2; % --> not yet initialized!!
        elseif P.j==P.jstar
            P.h1(:,P.j)         = P.h1star(:,P.j);
        end
        P.h1star(:,P.j)         = P.h1(:,P.j);
%% calcolo valori potenziale osmotico per il tempo di simulazione (ERROR)
        if W.iosm==1 && V.ifs>3
            P.IEC               = 0;        % boolean: store if ktec incremented
            %# ERROR: Antonio checks!! --> P.time(P.j)>EC.t(P.ktec)
            if P.time(P.j)>EC.t(P.ktec+1)
                P.ktec          = P.ktec+1;
                P.IEC           = 1;
            end
            P.ECstar(:,P.j)     = ( P.time(P.j)-EC.t(P.ktec) )      / ...
                              ( EC.t(P.ktec+1)-EC.t(P.ktec) )       * ...
                              ( EC.data(:,P.ktec+1)-EC.data(:,P.ktec) ) + ...
                                EC.data(:,P.ktec) ;
        end
%% attribuzione param ritenz e conducib ai vari nodi nei diversi strati
        for i=1:W.nz
            % curva ritenzione
            P.teta(i,P.j)       = fnteta( P.h1(i,P.j), P.sh, i );

            if or(P.sh.ifc==1,P.sh.ifc==3)
                P.kond(i,P.j)   = fncond( P.teta(i,P.j),P.sh, i );
            else
                P.kond(i,P.j)   = fncond( P.h1(i,P.j),  P.sh, i );
            end
            P.cap(i,P.j)        = fncap(  P.h1(i,P.j),  P.sh, i );

            if (and(W.iveg==1,W.itopvar==1))
                if P.Droot>0
                    P.dpt       = P.nodes.z(i);
                    if W.iosm==1 && V.ifs>3
                        P.op    = -P.ECstar(i,P.j)*360;
                    else
                        P.op    = 0;
                    end
                    if P.nodes.z(i)<P.Droot                    
                        P.sink(i,P.j)   = fnsink( P.h1(i,P.j), P, W, V );
                    else
                        P.sink(i,P.j)   = 0;
                    end
                else
                    P.sink(i,P.j)       = 0;
                end
            end        
        end
%% controllo sui flussi al contorno superiore
        if W.itopvar==1 && P.L==0 % L ==> tutto ok nel j precedente
            P.tq                    = B.top.thqstar(P.kk+1);

            if W.itbc==0 % flux
                W.qsurf             = B.top.hqstar(P.kk);
            elseif and(W.itbc==1,P.rnf==0) % potential
                W.hsurf             = B.top.hqstar(P.kk);
            end

            if W.ibotvar==1
                if W.ibbc==0
                    W.qbot          = B.bot.hqstar(P.kk);
                elseif W.ibbc==1
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
%% controllo tempi top boundary
        if W.itopvar==1
            if abs(P.time(P.j)-P.tq)<0.00001
                P.T                 = 1;
                P.L                 = 0;
            elseif P.time(P.j)>P.tq
                P.time(P.j)         = P.tq;
                W.dt              = P.tq-P.time(P.j-1);
                P.T                 = 1;
                P.L                 = 0;
            else
                P.T                 = 0;
                P.L                 = 1;
            end
        end
%% controllo sui tempi di stampa
        if abs(P.time(P.j)-W.tp(P.pp))<0.0001
            P.TT                    = 1;
        elseif P.time(P.j)>W.tp(P.pp)
            P.time(P.j)             = W.tp(P.pp);
            W.dt                  = W.tp(P.pp)-P.time(P.j-1);
            P.TT                    = 1;
        else
            P.TT                    = 0;
        end
%% flussi in superficie
        P.km_max(1,P.j)             = (P.kond(1,P.j)+P.sh.k0(1))/2;
        P.fluxsurf_max(P.j)         = -P.km_max(1,P.j)*((W.hsurfmax-P.h1(1,P.j))/(P.nodes.dz(1)/2)+1);% Darcy

        if W.itbc==0
            if and(and(and(W.itopvar==1,W.qsurf==0),W.iveg==1),W.itbc==0)
                W.hsurf             = 13.3*10^5*log(W.vpr);
                P.teta_hsurf        = fnteta( W.hsurf, P.sh, 1 );

                if or(P.sh.ifc==1,P.sh.ifc==3)
                    P.km(1,P.j)     = (P.kond(1,P.j) + fncond(P.teta_hsurf,P.sh,1))/2;
                else
                    P.km(1,P.j)     = (P.kond(1,P.j) + fncond(W.hsurf,P.sh,1))/2;
                end

                P.Emax              = -P.km(1,P.j) * ...
                                     ((W.hsurf-P.h1(1,P.j))/(P.nodes.dz(1)/2)+1);

                if P.Ep<P.Emax
                    W.qsurf         = P.Ep;
                    W.hsurf         = (3*P.h1(1,P.j)-P.h1(2,P.j))/2;
                    O.fluxsurf(1,P.j,mm) = W.qsurf;
                else
                    W.qsurf         = P.Emax;
                    O.fluxsurf(1,P.j,mm) = W.qsurf;
                end

            elseif and(W.itbc==0,W.qsurf<0)        
                W.hsurf             = (3*P.h1(1,P.j)-P.h1(2,P.j))/2;
                O.fluxsurf(1,P.j,mm)     = W.qsurf;
            elseif and(and(W.itbc==0,W.qsurf==0),W.iveg==0)
                W.hsurf             = (3*P.h1(1,P.j)-P.h1(2,P.j))/2;
                O.fluxsurf(1,P.j,mm)     = W.qsurf;
            end

            % Dall'estrapolazione lineare dei valori di h(1,P.j) e di
            % h(2,P.j) potrebbe risultare un W.hsurf>0 potendo tuttavia
            % ancora essere abs(W.qsurf)<abs(P.fluxsurf_max).
            % In questo caso P.rnf viene posto =1 % anche se alla fine
            % nella matrice RN compariranno valori positivi di runoff, che
            % ovviamente devono invece essere negativi.
            % Inoltre, W.itbc viene posto =1 in maniera da imporre un
            % potenziale al contorno superiore.
            if W.qsurf<0 && or( abs(W.qsurf)>=abs(P.fluxsurf_max(P.j)), ...
                                W.hsurf>=W.hsurfmax )
                P.rnf               = 1;
                W.itbc              = 1;
            end
        end

        if W.itbc==1
            if P.rnf==1
                W.hsurf             = W.hsurfmax;
                P.km_max(1,P.j)     = (P.kond(1,P.j)+P.sh.k0(1))/2;
                P.fluxsurf_max(P.j) = -P.km_max(1,P.j) * ...
                                  ((W.hsurfmax-P.h1(1,P.j))/(P.nodes.dz(1)/2)+1);
                O.fluxsurf(1,P.j,mm)= P.fluxsurf_max(P.j);
            elseif P.rnf==0
                P.teta_hsurf        = fnteta( W.hsurf, P.sh, 1 );
                if or(P.sh.ifc==1,P.sh.ifc==3)
                    P.km(1,P.j)     = ( P.kond(1,P.j) + ...
                                        fncond(P.teta_hsurf,P.sh,1) )/2;
                else
                    P.km(1,P.j)     = ( P.kond(1,P.j) + ...
                                        fncond(W.hsurf,P.sh,1) )/2;
                end
                O.fluxsurf(1,P.j,mm)     = -P.km(1,P.j) * ...
                                     ((W.hsurf-P.h1(1,P.j))/(P.nodes.dz(1)/2)+1);
            end
        end
%% flussi al fondo
        if W.ibbc==2
            W.hbot                  = P.h1(W.nz,P.j)-(W.grad-1)*P.nodes.dz(W.nz+1);
        end

        if W.ibbc==0
            O.fluxbot(1,P.j,mm)          = W.qbot;
            W.hbot(P.j)             = (P.h1(W.nz,P.j)-P.h1(W.nz-1,P.j))*...
                                       P.nodes.dz(W.nz+1)/P.nodes.dz(W.nz)+P.h1(W.nz,P.j);
        elseif or(W.ibbc==1,W.ibbc==2)
            P.teta_hbot             = fnteta( W.hbot, P.sh, W.nz);
            if or(P.sh.ifc==1,P.sh.ifc==3)
                P.kp(W.nz,P.j)      = ( P.kond(W.nz,P.j) + ...
                                        fncond(P.teta_hbot,P.sh,W.nz) )/2;
            else
                P.kp(W.nz,P.j)      = ( P.kond(W.nz,P.j) + ...
                                        fncond(W.hbot,P.sh,W.nz) )/2;
            end
            O.fluxbot(1,P.j,mm)          = -P.kp(W.nz,P.j) * ...
                                      ((P.h1(W.nz,P.j)-W.hbot)/P.nodes.dz(W.nz+1)+1);
        end
%% flussi ai nodi intermedi
        % Check errors for dz(?):
        %   > at P.flux(1,.)        --> dz(1)           *error
        %   > at P.flux(W.nz,.)     --> dz(end)         *error
        %   > at P.flux(2:W.nz-1,.) --> dz(2:W.nz-1)    *good!
        
        P.flux(1,P.j)               = -P.kond(1,P.j) * ...
                                   ((W.hsurf-P.h1(2,P.j))/(1.5*P.nodes.dz(1))+1);
        P.flux(W.nz,P.j)            = -P.kond(W.nz,P.j) * ...
                           ((P.h1(W.nz-1,P.j)-W.hbot)/(P.nodes.dz(W.nz+1)+P.nodes.dz(W.nz))+1);
        
        P.flux(2:W.nz-1,P.j)        = -P.kond(2:W.nz-1,P.j)     .*( ...
                                      ( P.h1((2:W.nz-1)-1,P.j)  -  ...
                                      P.h1((2:W.nz-1)+1,P.j) )  ./ ...
                                      (2*P.nodes.dz(2:W.nz-1))+1      );
%% calcolo flussi all'interfaccia
        P.flux(P.nodes.cumsum(2:W.nlay),P.j) = 2/3 * P.flux(P.nodes.cumsum(2:W.nlay)-1,P.j) +...
                                               1/3 * P.flux(P.nodes.cumsum(2:W.nlay)+2,P.j);
%% flussi cumulati in superficie ==> not used
%         if P.j==1
%             fluxsurfcum(P.j)        = O.fluxsurf(1,P.j,mm)*W.dt;
%         else
%             fluxsurfcum(P.j)        = O.fluxsurf(1,P.j,mm)*W.dt + ...
%                                       fluxsurfcum(P.j-1);
%         end
%% flussi cumulati ai nodi intermedi ==> not used
%         for i=1:W.nz
%             if P.j==1
%                 fluxcum(i,1)        = P.flux(i,1)*W.dt;
%             else
%                 fluxcum(i,P.j)      = P.flux(i,P.j)*W.dt + ...
%                                       fluxcum(i,P.j-1);
%             end
%         end
%% soluzione del sistema tridiagonale (see WARNING)
       P.h2                        = fnsyst( P, W );
%% calcolo capacita' metodo implicito (vedi Karvonen)
        for i=1:W.nz
            if abs(P.h2(i)-P.h1(i,P.j))>W.tolle1
                P.cap(i,P.j)       = ( fnteta( P.h2(i),P.sh,i) -  ...
                                        fnteta(P.h1(i,P.j),P.sh,i) ...
                                      ) / ( P.h2(i)-P.h1(i,P.j) ) ;
            end
        end
        % fnsyst con il nuovo P.cap:
        O.h22(:,P.j,mm)            = fnsyst( P, W );
%% calcolo capacita' al tempo medio
        % W.tolle1=tolleranza per la convergenza della capacit�
        P.SS                       = 0;
        P.niter                    = 0; % Nj per convergenza capacit�.
        while and(max(abs(O.h22(:,P.j,mm)-P.h2))>W.tolle1,P.SS==0)
            for i=1:W.nz            
                if abs(O.h22(i,P.j,mm)-P.h2(i))>W.tolle1        
                    P.cap(i,P.j)    = ( fnteta(O.h22(i,P.j,mm),P.sh,i) -   ...
                                        fnteta(P.h1(i,P.j),P.sh,i)     ...
                                      ) / ( O.h22(i,P.j,mm)-P.h1(i,P.j) )  ;
                end
            end
            P.h2                    = O.h22(:,P.j,mm);
            O.h22(:,P.j,mm)         = fnsyst( P, W );

            P.niter                = P.niter+1;
            % Migliorare soluzione, o comunque salvare in M.nnc l'eventuale simulazione saltata
            if P.niter>10
                P.SS                = 1;
                if W.dt==W.dtmin
                    P.LL            = 0;
                elseif W.dt>W.dtmin
                    P.LL            = 1;
                    W.dt          = W.dt/3;
                    P.time(P.j)     = P.time(P.j)-2*W.dt;
                    if W.dt<=W.dtmin
                        P.time(P.j) = P.time(P.j) - W.dt + W.dtmin;
                        W.dt      = W.dtmin;
                    end
                    % Ripristino flussi o potenziali nel caso di
                    % annullamento dell'iterazione per riduzione del
                    % W.dt.
                    % Serve solo nel caso di P.LL=1, quando il ripristino
                    % non pu� avvenire con le righe uguali poste alla fine
                    % del blocco per la simulazione del trasporto di
                    % soluti.
                    if and(W.itopvar==1,W.itbc==0)
                        W.qsurf     = B.top.hqstar(P.kk);
                    elseif and(and(W.itopvar==1,W.itbc==1),P.rnf==0)
                        W.hsurf     = B.top.hqstar(P.kk);
                    elseif and(W.itbc==1,P.rnf==1)
                        W.qsurf     = B.top.hqstar(P.kk);
                    end
                    % Ripristino tempo di lettura delle EC.data nel caso in
                    % cui l'annullamento della simulazione avvenga a
                    % cavallo del cambio di EC.t
                    if W.iosm==1 && V.ifs>3 && P.IEC==1
                        P.ktec      = P.ktec-1;
                    end
                end
            end
        end
        
        if P.niter<=3
            W.dt                  = W.multmin * W.dt;
            if W.dt>W.dtmax
                W.dt              = W.dtmax;
            end
            P.LL                  = 0;
            
        elseif and(P.niter>=4,P.niter<=10)
            W.dt                  = W.multmax * W.dt;
            if W.dt<W.dtmin
                W.dt              = W.dtmin;
            end
            P.LL                  = 0;
        end
%% ?? define cell ??
        if P.LL==0
%% trasporto soluti CDE
            if W.isol==2
                P = solute_transport_ADE_N( P, W, S, B, O, mm );
            end
%% restore qsurf & hsurf -- check with Antonio!!
% Serve a ripristinare il valore di W.qsurf che potrebbe essere stato
% cambiato nella routine per l'evaporazione W.qsurf=P.Ep oppure
% W.qsurf=P.Emax.
% Lo stesso vale per W.hsurf che potrebbe essere stato cambiato nella
% routine per W.itbc=1 W.hsurf=(3*P.h1(1,P.j)-P.h1(2,P.j))/2 oppure
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
    end% P.time(P.j)<W.tmax
%% MONTECARLO --END--
end% mm=1:M.nvp
%% SAVE -- incomplete
if W.MTCL == 0
    multilayer_save( O )
elseif W.MTCL == 1
    multilayer_save_mcs.m
end