%% initialization
P.A                             = NaN(P.Nj,W.nz+2);
P.R                             = NaN(P.Nj,5);
P.Q                             = NaN(P.Nj,4);
P.time                          = NaN(P.Nj,1);
P.h1                            = NaN(W.nz,2,P.Nj); % **
P.h1star                        = NaN(W.nz,2,P.Nj); % **
P.j                             = NaN;
P.jstar                         = NaN;
P.C1                            = NaN(W.nz,2,P.Nj);
P.C1star                        = NaN(W.nz,2,P.Nj);
P.C2                            = NaN(W.nz,2,P.Nj);
P.S1                            = NaN(W.nz,2,P.Nj);
P.S2                            = NaN(W.nz,2,P.Nj);
P.ECstar                        = NaN(W.nz,P.Nj);
P.istar                         = NaN(W.nlay+1,1);
P.dz                            = NaN(W.nz,1);
P.dap                           = NaN(W.nz,1);
P.tetas                         = NaN(W.nz,1);
P.tetar                         = NaN(W.nz,1);
P.alfrs                         = NaN(W.nz,1);
P.fi                            = NaN(W.nz,1);
P.alfvg                         = NaN(W.nz,1);
P.en                            = NaN(W.nz,1);
P.alfvg2                        = NaN(W.nz,1);
P.en2                           = NaN(W.nz,1);
P.ifr                           = NaN(W.nz,1);
P.k0                            = NaN(W.nz,1);
P.k0macr                        = NaN(W.nz,1);
P.bita                          = NaN(W.nz,1);
P.bita2                         = NaN(W.nz,1);
P.ifc                           = NaN(W.nz,1);
P.teta                          = NaN(W.nz,P.Nj);
P.tetafc                        = NaN(W.nz,1);
P.kond                          = NaN(W.nz,P.Nj);
P.cap                           = NaN(W.nz,P.Nj);
P.dpt                           = NaN;
P.op                            = NaN;
P.sink                          = NaN(W.nz,P.Nj);
P.Cinput                        = NaN(2,1);
P.km_max                        = NaN(1,P.Nj);
P.fluxsurf_max                  = NaN(1,P.Nj);
P.teta_hsurf                    = NaN;
P.km                            = NaN(1,P.Nj);
P.Emax                          = NaN;
P.fluxsurf                      = NaN(1,P.Nj);
P.runoff                        = NaN(1,P.Nj);
P.kp                            = NaN(W.nz,P.Nj);
P.fluxbot                       = NaN(1,P.Nj);
P.teta_hbot                     = NaN;
P.flux                          = NaN(W.nz,P.Nj);
P.k                             = NaN;
P.h2                            = NaN(W.nz,1);
P.h22                           = NaN(W.nz,1);
P.SS                            = NaN;
P.niter                         = NaN;
%% MONTECARLO --START--
for mm=1:M.nvp
    progress_bar(mm,M.nvp,'multilayer prog')
%% Assign vars at soil layers in case of stochastic simulation
    if W.MTCL==1
        % Assign the stochastic values to the variables to be simulated
        % using Montecarlo:
        for nvars=1:length(M.list)
            % ANTONIO: IS THAT FINE? IS WHAT WE NEED?
            eval( [M.list{ nvars } ' = M.data(M.combinations(',num2str(mm),',:),',num2str(nvars),')'';'] )
        end
    end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% fine Routine MONTECARLO MODULE I 
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%% Definition of water movement field
% Definizione del campo di moto
    % definire meglio i seguenti(!!!):
    %   > P.istar:    numero del nodo immediatamente sopra l'interfaccia
    %   > P.istar(1): nodo alla superficie
    % -----perch� il seguente pezzo non st� in _conf.m??------
    P.istar(1)          = 0;
    P.istar(W.nlay+1)   = W.nz;
    W.zint              = [0,W.zint]; % perch� non va inserito lo zero in config?
    P.zmax              = W.zint(W.nlay+1);
    P.d_z               = NaN(1,W.nlay);
    P.z                 = NaN(1,P.istar(end)+1);
    % --------------------------------------------------------

    % Multistrato:
    if W.nlay>1
        for i=2:W.nlay
            % P.istar: nodo all'interfaccia.
            P.istar(i)   = P.istar(i-1) + int32( W.nz*(W.zint(i)-W.zint(i-1))/P.zmax );
            P.d_z(i-1)   = (W.zint(i)-W.zint(i-1)) / (P.istar(i)-P.istar(i-1));
        end
        % calcola lo spessore dell'internodo di ciascuno strato:
        P.d_z(W.nlay)   = (W.zint(W.nlay+1)-W.zint(W.nlay))/(P.istar(W.nlay+1)-P.istar(W.nlay));
        P.z(1)          = P.d_z(1)/2;
        P.dztop         = P.z(1);
        for j=2:W.nlay+1
           for i=P.istar(j-1)+2:P.istar(j)+1
               P.z(i)   = P.z(i-1)+P.d_z(j-1);
           end
        end

    % Monostrato:
    else
        P.d_z(1)        = P.zmax/W.nz;
        P.z(1)          = P.d_z(1)/2; % P.z ha un elemento in meno rispetto a Multistrato! E' corretto?? 
        P.dztop         = P.z(1);
        for i=2:W.nz
            P.z(i)      = P.z(i-1)+P.d_z(1);
        end
    end

    % Definisce lo spessore dell'ultimo dz:
    P.dzbot             = 2*(P.zmax-P.z(W.nz));
%% Create printing matrices -- put in initialization --
    
    W.dt                = W.dtin;
    
    % Perch� duplicare A spasmodicamente?
    % Ad esempio, perch� B esiste se � uguale ad A?? --> vedi anche dopo...
    P.A(1:2,:)          = [ -1:1:W.nz; [-1,0,P.z(1:W.nz)] ];
%     B                   = A;
%     C                   = A;
%     DNH                 = A;
%     DNO                 = A;
%     E                   = A;
%     F                   = A;
%     G                   = A;
    P.R(1,:)            = 1:5;
%     RN                  = R;
    P.Q(1,:)            = 1:4;
%     Qtb                 = Q;

    % inizializzazione contatori e tempo di simulazione
    P.j                 = 1; % mettiamo i contatori in struct array specifica come "C"
    P.time(1)           = W.dt;
    P.kk                = 1;
    P.pp                = 1;
    P.L                 = 0;
    P.LL                = 0;
%     P.RR                = 1; % never used
%     P.FF                = 0; % never used
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
%% definition of retention & conductivity pars at nodes
    for k=1:W.nlay
        for i=P.istar(k)+1:P.istar(k+1)
            P.dz(i)             = P.d_z(k);
            P.dap(i)            = W.d_ap(k);
            P.tetas(i)          = W.tetas(k);
            P.tetar(i)          = W.tetar(k);
            P.alfrs(i)          = W.alfrs(k);
            P.fi(i)             = W.fi(k);
            P.alfvg(i)          = W.alfvg(k);
            P.en(i)             = W.en(k);
            P.alfvg2(i)         = W.alfvg2(k);
            P.en2(i)            = W.en2(k);
            P.ifr(i)            = W.ifr(k);
            P.k0(i)             = W.k0(k);
            P.k0macr(i)         = W.k0macr(k);
            P.bita(i)           = W.bita(k);
            P.bita2(i)          = W.bita2(k);
            P.ifc(i)            = W.ifc(k);

            % Calcolo del teta alla field capacity da utilizzare nel
            % modulo soluti per il calcolo del fattore di riduzione del
            % coefficiente di mineralizzazione.
            P.tetafc(i)         = fnteta( W.hfc, P, i );
        end
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
%% load initial conc. & update conc. for CDE transport model

    % inserito nel modulo suo!!
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
            P.teta(i,P.j)       = fnteta( P.h1(i,P.j), P, i );

            if or(P.ifc==1,P.ifc==3)
                P.kond(i,P.j)   = fncond( P.teta(i,P.j), P, i );
            else
                P.kond(i,P.j)   = fncond( P.h1(i,P.j), P, i );
            end
            P.cap(i,P.j)        = fncap(  P.h1(i,P.j), P, i );

            if (and(W.iveg==1,W.itopvar==1))
                if P.Droot>0
                    P.dpt       = P.z(i);
                    if W.iosm==1 && V.ifs>3
                        P.op    = -P.ECstar(i,P.j)*360;
                    else
                        P.op    = 0;
                    end
                    if P.z(i)<P.Droot                    
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
        if W.itopvar==1 && P.L==0
            P.tq                    = B.top.thqstar(P.kk+1);

            if W.itbc==0
                W.qsurf             = B.top.hqstar(P.kk);
            elseif and(W.itbc==1,P.rnf==0)
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
        P.km_max(1,P.j)             = (P.kond(1,P.j)+P.k0(1))/2;
        P.fluxsurf_max(P.j)         = -P.km_max(1,P.j)*((W.hsurfmax-P.h1(1,P.j))/(P.dz(1)/2)+1);

        if W.itbc==0
            if and(and(and(W.itopvar==1,W.qsurf==0),W.iveg==1),W.itbc==0)
                W.hsurf             = 13.3*10^5*log(W.vpr);
                P.teta_hsurf        = fnteta( W.hsurf, P, 1 );

                if or(P.ifc==1,P.ifc==3)
                    P.km(1,P.j)     = (P.kond(1,P.j) + fncond(P.teta_hsurf,P,1))/2;
                else
                    P.km(1,P.j)     = (P.kond(1,P.j) + fncond(W.hsurf,P,1))/2;
                end

                P.Emax              = -P.km(1,P.j) * ...
                                     ((W.hsurf-P.h1(1,P.j))/(P.dz(1)/2)+1);

                if P.Ep<P.Emax
                    W.qsurf         = P.Ep;
                    W.hsurf         = (3*P.h1(1,P.j)-P.h1(2,P.j))/2;
                    P.fluxsurf(P.j) = W.qsurf;
                else
                    W.qsurf         = P.Emax;
                    P.fluxsurf(P.j) = W.qsurf;
                end

            elseif and(W.itbc==0,W.qsurf<0)        
                W.hsurf             = (3*P.h1(1,P.j)-P.h1(2,P.j))/2;
                P.fluxsurf(P.j)     = W.qsurf;
            elseif and(and(W.itbc==0,W.qsurf==0),W.iveg==0)
                W.hsurf             = (3*P.h1(1,P.j)-P.h1(2,P.j))/2;
                P.fluxsurf(P.j)     = W.qsurf;
            end

            % Dall'estrapolazione lineare dei valori di h(1,P.j) e di
            % h(2,P.j) potrebbe risultare un W.hsurf>0 potendo tuttavia
            % ancora essere abs(W.qsurf)<abs(P.fluxsurf-max).
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
                P.km_max(1,P.j)     = (P.kond(1,P.j)+P.k0(1))/2;
                P.fluxsurf_max(P.j) = -P.km_max(1,P.j) * ...
                                  ((W.hsurfmax-P.h1(1,P.j))/(P.dz(1)/2)+1);
                P.fluxsurf(P.j)     = P.fluxsurf_max(P.j);
            elseif P.rnf==0
                P.teta_hsurf        = fnteta( W.hsurf, P, 1 );
                if or(P.ifc==1,P.ifc==3)
                    P.km(1,P.j)     = ( P.kond(1,P.j) + ...
                                        fncond(P.teta_hsurf,P,1) )/2;
                else
                    P.km(1,P.j)     = ( P.kond(1,P.j) + ...
                                        fncond(W.hsurf,P,1) )/2;
                end
                P.fluxsurf(P.j)     = -P.km(1,P.j) * ...
                                     ((W.hsurf-P.h1(1,P.j))/(P.dz(1)/2)+1);
            end
        end
%% flussi al fondo
        if W.ibbc==2
            W.hbot                  = P.h1(W.nz,P.j)-(W.grad-1)*P.dzbot;
        end

        if W.ibbc==0
            P.fluxbot(P.j)          = W.qbot;
            W.hbot(P.j)             = (P.h1(W.nz,P.j)-P.h1(W.nz-1,P.j))*...
                                       P.dzbot/P.dz(W.nz)+P.h1(W.nz,P.j);
        elseif or(W.ibbc==1,W.ibbc==2)
            P.teta_hbot             = fnteta( W.hbot, P, W.nz);
            if or(P.ifc==1,P.ifc==3)
                P.kp(W.nz,P.j)      = ( P.kond(W.nz,P.j) + ...
                                        fncond(P.teta_hbot,P,W.nz) )/2;
            else
                P.kp(W.nz,P.j)      = ( P.kond(W.nz,P.j) + ...
                                        fncond(W.hbot,P,W.nz) )/2;
            end
            P.fluxbot(P.j)          = -P.kp(W.nz,P.j) * ...
                                      ((P.h1(W.nz,P.j)-W.hbot)/P.dzbot+1);
        end
%% flussi ai nodi intermedi
        % Check errors for dz(?):
        %   > at P.flux(1,.)        --> dz(1)           *error
        %   > at P.flux(W.nz,.)     --> dz(end)         *error
        %   > at P.flux(2:W.nz-1,.) --> dz(2:W.nz-1)    *good!
        
        P.flux(1,P.j)               = -P.kond(1,P.j) * ...
                                   ((W.hsurf-P.h1(2,P.j))/(1.5*P.dz(1))+1);
        P.flux(W.nz,P.j)            = -P.kond(W.nz,P.j) * ...
                           ((P.h1(W.nz-1,P.j)-W.hbot)/(P.dzbot+P.dz(W.nz))+1);
        
        P.flux(2:W.nz-1,P.j)        = -P.kond(2:W.nz-1,P.j)     .*( ...
                                      ( P.h1((2:W.nz-1)-1,P.j)  -  ...
                                      P.h1((2:W.nz-1)+1,P.j) )  ./ ...
                                      (2*P.dz(2:W.nz-1))+1      );
%% calcolo flussi all'interfaccia
        P.k                         = 2:W.nlay;
        P.flux(P.istar(P.k),P.j)    = 2/3 * P.flux(P.istar(P.k)-1,P.j) +...
                                      1/3 * P.flux(P.istar(P.k)+2,P.j);
%% flussi cumulati in superficie ==> not used
%         if P.j==1
%             fluxsurfcum(P.j)        = P.fluxsurf(P.j)*W.dt;
%         else
%             fluxsurfcum(P.j)        = P.fluxsurf(P.j)*W.dt + ...
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
       P.h2                         = fnsyst( P, W );
%% calcolo capacita' metodo implicito (vedi Karvonen)
        for i=1:W.nz
            if abs(P.h2(i)-P.h1(i,P.j))>W.tolle1
                P.cap(i,P.j)        = ( fnteta( P.h2(i),P,i) -  ...
                                        fnteta(P.h1(i,P.j),P,i) ...
                                      ) / ( P.h2(i)-P.h1(i,P.j) ) ;
            end
        end
        % fnsyst con il nuovo P.cap:
        P.h22                       = fnsyst( P, W );
%% calcolo capacita' al tempo medio
        % W.tolle1=tolleranza per la convergenza della capacit�
        P.SS                        = 0;
        P.niter                    = 0; % Nj per convergenza capacit�.
        while and(max(abs(P.h22-P.h2))>W.tolle1,P.SS==0)
            for i=1:W.nz            
                if abs(P.h22(i)-P.h2(i))>W.tolle1        
                    P.cap(i,P.j)    = ( fnteta(P.h22(i),   P,i) -   ...
                                        fnteta(P.h1(i,P.j),P,i)     ...
                                      ) / ( P.h22(i)-P.h1(i,P.j) )  ;
                end
            end
            P.h2                    = P.h22;
            P.h22                   = fnsyst( P, W );

            P.niter                = P.niter+1;
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
                P = solute_transport_ADE_N( P, W, S, B );
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
%             if and(W.itopvar==1,W.itbc==0)
%                 W.qsurf=B.top.hqstar(P.kk);
%             elseif and(and(W.itopvar==1,W.itbc==1),P.rnf==0)
%                 W.hsurf=B.top.hqstar(P.kk);
%             elseif and(W.itbc==1,P.rnf==1)
%                 W.qsurf=B.top.hqstar(P.kk);
%             end
            W.qsurf     = B.top.hqstar(P.kk);
            W.hsurf     = B.top.hqstar(P.kk);
%% print potentials -- check with Antonio!!

            % NOTA:
            % Per quello che ho capito le matrici importanti sono le
            % seguenti:
            %   > P.C2          --> concentrazioni soluti [nodi x tempi x soluti]
            %   > P.h22         --> flussi ai nodi intermedi ?? [nodi x tempi] 
            %   > P.fluxsurf    --> flusso al contorno superiore ??;
            %   > P.fluxbot     --> flusso al contorno inferiore ??;
            %   > W.qsurf       --> flusso imposto al contorno superiore ??;
            %   > P.runoff      --> runoff;
            %   > in una implementazione futura piu' integrata anche il
            %     runon
            %
            % I farei in questo modo: conserverei solo le matrici utili
            % elencate sopra, tal quali e senza rimaneggiamenti. inoltre in
            % Initialization prevederei una preallocazione considerando
            % anche il numero di simulazioni montecarlo cos da mantenere
            % tutto quanto - ridotto ai minimi termini - in memoria RAM, e
            % solo alla fine di tutto il run si salva quello che il
            % programma DEVE salvare (le info da salvare possono anche
            % essere configurate dall'utente, per cui si prevede una
            % sezione in 'conf' apposita).

            % stampa tutti i potenziali nella matrice P.A
            P.A(P.j+2,:)    = [P.j;P.time(P.j);P.h22]';
            P.runoff(P.j)     = W.qsurf-P.fluxsurf(P.j);
            P.R(P.j+1,:)    = [P.j,P.time(P.j),W.qsurf,P.fluxsurf(P.j),P.runoff(P.j)];
            P.Q(P.j+1,:)    = [P.j,P.time(P.j),P.fluxsurf(P.j),P.fluxbot(P.j)];

            if P.rnf==1
                RN=[RN;P.R(P.j+1,:)];
            end

            % stampa le variabili per i soli tempi di stampa e aggiorna il
            % contatore per il tempo di stampa
            if P.TT==1
                % lascialo!! %
                P.pp                        = P.pp+1;
                % lascialo!! %
                
%                 B=[B;P.A(P.j+2,:)]; % perch� B esiste se � uguale ad A??
%                 Qtb=[Qtb;P.Q(P.j+1,:)];
                for i=1:W.nz
                    % a cosa serve il tetaout rispetto ai teta gia'
                    % calcolati?
                    tetaout(i)              = fnteta(P.h22(i),P,i);

                    if and(W.iveg==1,W.itopvar==1)
                        if P.Droot>0
                            P.dpt           = P.z(i);
                            if W.iosm==1
                                if V.ifs>3
                                    P.op    = -P.ECstar(i,P.j)*360;
                                end
                            else
                                P.op        = 0;
                            end
                            if P.z(i)<P.Droot                            
                                sinkout(i)  = fnsink( P.h1(i,P.j), P, W, V );
                            else
                                sinkout(i)  = 0;
                            end
                        else
                            sinkout(i)      = 0;
                        end
                    end
                    capout(i)               = fncap(  P.h1(i,P.j), P, i );
                    fluxout(i)              = P.flux(i,P.j);
                end
                C=[C;P.j,P.time(P.j),tetaout];
                if or(W.isol==1,W.isol==2)
                    DNH=[DNH;P.j,P.time(P.j),C2out_NH];
                    DNO=[DNO;P.j,P.time(P.j),C2out_NO];                   
                end
                E=[E;P.j,P.time(P.j),sinkout];
                F=[F;P.j,P.time(P.j),capout];
                G=[G;P.j,P.time(P.j),fluxout];
            end
%% controllo di W.itbc -- controllare con Antonio!!
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
            if P.T==1               % che vuol dire?
                P.kk=P.kk+1;                
                 if W.itopvar==1    % che vuol dire?
                     if and(P.rnf==1,abs(B.top.hqstar(P.kk))<abs(B.top.hqstar(P.kk-1))) % che vuol dire?
                        P.rnf=0;   
                        W.itbc=0;
                     end
                 end
            end
%% update time of simulation [ok!]
            P.j         = P.j+1;
            P.time(P.j) = P.time(P.j-1)+W.dt;
        end% if P.LL=0    
    end% P.time(P.j)<W.tmax
%% Traspone matrice -- DELETE!!!
    P.A=P.A';
    B=B';
    C=C';
    DNH=DNH';
    DNO=DNO';
    E=E';
    F=F';
    G=G';
    RN=RN';
    Qtb=Qtb';
%% da scommentare in caso di simulazione MONTECARLO
    if W.MTCL==1

        %%%%MONTECARLO TETA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if P.LL==0
            [r,c] = size(C);

            %OPEN file
            fid = fopen( strcat(proj.path,'teta_print',num2str(mm),'.dat'), 'w' );
            %cicli loop
            for riga=1:r
                for colonna = 1:c
                    v = C(riga,colonna);
                    if  colonna == c
                        fprintf(fid,'%f\n',v);
                    else
                        fprintf(fid,'%f ',v);
                    end
                end
            end
            state = fclose(fid);
            close;

            M.tetasum=C+M.tetasum;
            M.tetasumSQ=(C.^2)+M.tetasumSQ;


        %%%%MONTECARLO CONCENTRAZIONI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [rr,cc] = size(D);

            %OPEN file
            fid = fopen( strcat(proj.path,'conc_print',num2str(mm),'.dat'), 'ww' );
            %cicli loop
            for riga=1:rr
                for colonna = 1:cc
                    vv = D(riga,colonna);
                    if  colonna == cc
                        fprintf(fid,'%f\n',vv);
                    else
                        fprintf(fid,'%f ',vv);
                    end
                end
            end
            state = fclose(fid);
            close;

            M.concsum=D+M.concsum;
            M.concsumSQ=(D.^2)+M.concsumSQ;

        %%%%MONTECARLO FLUSSI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [rrr,ccc] = size(G);

            %OPEN file
            fid = fopen([strcat(proj.path,'flux_i_j',num2str(mm),'.dat')],'ww');
            %cicli loop
            for riga=1:rr
                for colonna = 1:ccc
                    vvv = G(riga,colonna);
                    if  colonna == ccc
                        fprintf(fid,'%f\n',vvv);
                    else
                        fprintf(fid,'%f ',vvv);
                    end
                end
            end
            state = fclose(fid);
            close;

            M.fluxsum=G+M.fluxsum;
            M.fluxsumSQ=(G.^2)+M.fluxsumSQ;

        %fine secondo if P.LL=0
        end

        if and(W.itopvar==1,W.itbc==1)
            W.itbc=0;
        end
    end
%% MONTECARLO --END--
end% mm=1:M.nvp
%% SAVE -- edit with Antonio

multilayer_save( P )