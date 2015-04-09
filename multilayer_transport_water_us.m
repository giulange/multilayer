% multilayer_transport_water_su
% 
% DESCRIPTION
%   This script simulates the transport of water in unsaturated/saturated
%   soils at a time step specified by the main "multilayer" routine.
%   It solves the Richard's equation at each compartment of the soil
%   system.
%   ...top boundary implementation
%   ...bottom boundary implementation
%   ...more hydrological details on what is performed...
%   It can handle a fixed and a variable compartment height.
%   ...other characteristics of the model implemented.

%% lettura h iniziali o aggiornamento h -- check at the end!
% lettura h iniziali e aggiornamento h
if P.j==1
    P.h1                    = W.hin;
elseif P.j>P.jstar  % if simulation runs normally
    P.h1                    = P.h2; % --> not yet initialized!!
elseif P.j==P.jstar % if simulation at previous step didn't go fine!
    P.h1                    = P.h1star;
end
P.h1star                    = P.h1;
%% calc. param. ritenz. e conducib. ai vari nodi nei diversi strati
for i=1:P.nz
    % Ritenzione:
    P.teta(i)               = fnteta( P.h1(i), P.sh, i );
    % Conducibilità:
    if P.sh.ifc(i)==1 || P.sh.ifc(i)==3
        P.kond(i)           = fncond( P.teta(i),P.sh, i );
    else
        P.kond(i)           = fncond( P.h1(i),  P.sh, i );
    end
    % Capacità:
    P.cap(i)                = fncap(  P.h1(i),  P.sh, i );
    % Sink:
    if W.iveg==1 && W.itopvar==1
        if P.Droot>0
            P.dpt           = P.nodes.z(i);
            if W.iosm==1 && V.ifs>3
                P.op        = -P.ECstar(i)*360;
            else
                P.op        = 0;
            end
            if P.nodes.z(i) < P.Droot                    
                P.sink(i)   = fnsink( P.h1(i), P, W, V );
            else
                P.sink(i)   = 0;
            end
        else
            P.sink(i)       = 0;
        end
    end
end
%% FLUX :: TOP boundary
P.km_max(P.j)               = (P.kond(1)+P.sh.k0(1))/2;
P.fluxsurf_max(P.j)         = -P.km_max(P.j)*((W.hsurfmax-P.h1(1))/(P.nodes.dz(1)/2)+1);% Darcy

switch W.itbc
    case 0
    if and(and(and(W.itopvar==1,W.qsurf==0),W.iveg==1),W.itbc==0) % delete W.itbc=0
        W.hsurf             = 13.3 * 10^5 * log(W.vpr);
        P.teta_hsurf        = fnteta( W.hsurf, P.sh, 1 );

        if or(P.sh.ifc==1,P.sh.ifc==3)
            P.km(P.j)       = (P.kond(1) + fncond(P.teta_hsurf,P.sh,1))/2;
        else
            P.km(P.j)       = (P.kond(1) + fncond(W.hsurf,P.sh,1))/2;
        end
        P.Emax              = -P.km(P.j) * ((W.hsurf-P.h1(1))/(P.nodes.dz(1)/2)+1);

        if P.Ep<P.Emax
            W.qsurf         = P.Ep;
            W.hsurf         = (3*P.h1(1)-P.h1(2))/2;
            O.fluxsurf(1,P.j,mm) = W.qsurf; % porta fuori (ma prima controlla che sia qsurf e non hsurf!!
        else
            W.qsurf         = P.Emax;
            O.fluxsurf(1,P.j,mm) = W.qsurf; % porta fuori
        end

    elseif and(W.itbc==0,W.qsurf<0) % delete W.itbc=0
        W.hsurf             = (3*P.h1(1)-P.h1(2))/2;
        O.fluxsurf(1,P.j,mm)= W.qsurf;
    elseif and(and(W.itbc==0,W.qsurf==0),W.iveg==0) % delete W.itbc=0
        W.hsurf             = (3*P.h1(1)-P.h1(2))/2;
        O.fluxsurf(1,P.j,mm)= W.qsurf;
    end

    % Dall'estrapolazione lineare dei valori di h(1,P.j) e di
    % h(2,P.j) potrebbe risultare un W.hsurf>0 potendo tuttavia
    % ancora essere abs(W.qsurf)<abs(P.fluxsurf_max).
    % In questo caso P.rnf viene posto =1 % anche se alla fine
    % nella matrice RN compariranno valori positivi di runoff, che
    % ovviamente devono invece essere negativi.
    % Inoltre, W.itbc viene posto =1 in maniera da imporre un
    % potenziale al contorno superiore.
    if W.qsurf<0 && or( abs(W.qsurf)>=abs(P.fluxsurf_max(P.j)),W.hsurf>=W.hsurfmax )
        P.rnf               = 1;
        W.itbc              = 1;
    end

    case 1
    if P.rnf==1
        W.hsurf             = W.hsurfmax;
        P.km_max(P.j)       = (P.kond(1)+P.sh.k0(1))/2;
        P.fluxsurf_max(P.j) = -P.km_max(P.j) * ...
                          ((W.hsurfmax-P.h1(1))/(P.nodes.dz(1)/2)+1);
        O.fluxsurf(1,P.j,mm)= P.fluxsurf_max(P.j);
    elseif P.rnf==0
        P.teta_hsurf        = fnteta( W.hsurf, P.sh, 1 );
        if or(P.sh.ifc==1,P.sh.ifc==3)
            P.km(P.j)       = ( P.kond(1) + ...
                                fncond(P.teta_hsurf,P.sh,1) )/2;
        else
            P.km(P.j)       = ( P.kond(1) + ...
                                fncond(W.hsurf,P.sh,1) )/2;
        end
        O.fluxsurf(1,P.j,mm)= -P.km(P.j) * ...
                             ((W.hsurf-P.h1(1))/(P.nodes.dz(1)/2)+1);
    end
    
    otherwise
    % nothing to do!
    
end
%% FLUX :: BOTTOM boundary
if W.ibbc==2
    W.hbot                  = P.h1(P.nz)-(W.grad-1)*P.nodes.dz(P.nz+1);
end

%         switch W.ibbc
%             case 0        
if W.ibbc==0
    O.fluxbot(1,P.j,mm)          = W.qbot;
    W.hbot(P.j)             = (P.h1(P.nz)-P.h1(P.nz-1))*...
                               P.nodes.dz(P.nz+1)/P.nodes.dz(P.nz)+P.h1(P.nz);
%             case {1,2}
elseif or(W.ibbc==1,W.ibbc==2)
    P.teta_hbot             = fnteta( W.hbot, P.sh, P.nz);
    if or(P.sh.ifc==1,P.sh.ifc==3)
        P.kp(P.nz)          = ( P.kond(P.nz) + fncond(P.teta_hbot,P.sh,P.nz) )/2;
    else
        P.kp(P.nz)          = ( P.kond(P.nz) + fncond(W.hbot,P.sh,P.nz) )/2;
    end
    O.fluxbot(1,P.j,mm)     = -P.kp(P.nz) * ((P.h1(P.nz)-W.hbot)/P.nodes.dz(P.nz+1)+1);
end
%% FLUX :: INTERMEDIATE nodes -- errors? check with Antonio
% Check errors for dz(?):
%   > at P.flux(1,.)        --> dz(1)           *error
%       - denominator: 1.5? +1? dz(1) soltanto?
%   > at P.flux(P.nz,.)     --> dz(end)         *error
%       - denominator: +1?
%   > at P.flux(2:P.nz-1,.) --> dz(2:P.nz-1)    *good!

P.flux(1)                   = -P.kond(1) * ((W.hsurf-P.h1(2)) /(1.5*P.nodes.dz(1))+1);
P.flux(P.nz)                = -P.kond(P.nz) * ((P.h1(P.nz-1)-W.hbot) ...
                              /(P.nodes.dz(P.nz+1)+P.nodes.dz(P.nz))+1);
P.flux(2:P.nz-1)            = -P.kond(2:P.nz-1).*( ...
                              ( P.h1((2:P.nz-1)-1)-P.h1((2:P.nz-1)+1) ) ...
                              ./ (2*P.nodes.dz(2:P.nz-1))+1 );
%% calcolo flussi all'interfaccia -- check with Antonio
% Definiamo l'obiettivo, poi valutiamo come implementare.
% P.k = 2:W.nlay;

% why +2 and not +1??
P.flux(P.nodes.cumsum(2:W.nlay)) = 2/3 * P.flux(P.nodes.cumsum(2:W.nlay)-1) +...
                                   1/3 * P.flux(P.nodes.cumsum(2:W.nlay)+2);
%% flussi cumulati in superficie ==> not used
%         if P.j==1
%             fluxsurfcum(P.j)        = O.fluxsurf(1,P.j,mm)*W.dt;
%         else
%             fluxsurfcum(P.j)        = O.fluxsurf(1,P.j,mm)*W.dt + ...
%                                       fluxsurfcum(P.j-1);
%         end
%% flussi cumulati ai nodi intermedi ==> not used -- error on P.j
%         for i=1:P.nz
%             if P.j==1
%                 fluxcum(i,1)        = P.flux(i,1)*W.dt;
%             else
%                 fluxcum(i,P.j)      = P.flux(i,P.j)*W.dt + ...
%                                       fluxcum(i,P.j-1);
%             end
%         end
%% soluzione del sistema tridiagonale (see WARNING) -- obiettivo?
% Obiettivo?
P.h2                        = fnsyst( P, W );
%% calcolo capacita' metodo implicito (vedi Karvonen) -- obietivo?
for i=1:P.nz
    if abs(P.h2(i)-P.h1(i))>W.tolle1
        P.cap(i)            = ( fnteta( P.h2(i),P.sh,i) - ...
                            fnteta(P.h1(i),P.sh,i) ) / ( P.h2(i)-P.h1(i) );
    end
end
% fnsyst con il nuovo P.cap:
O.h22(:,P.j,mm)             = fnsyst( P, W );
%         figure(13),hold on, plot(squeeze(O.h22(:,P.j,mm))), hold off;
%% calcolo capacita' al tempo medio -- obiettivo?
% W.tolle1=tolleranza per la convergenza della capacità
P.SS                        = 0;
P.niter                     = 0; % Nj per convergenza capacità.
while and(max(abs(O.h22(:,P.j,mm)-P.h2))>W.tolle1,P.SS==0)
    for i=1:P.nz            
        if abs(O.h22(i,P.j,mm)-P.h2(i))>W.tolle1        
            P.cap(i)        = ( fnteta(O.h22(i,P.j,mm),P.sh,i) - ...
                                fnteta(P.h1(i),P.sh,i) ) / ( O.h22(i,P.j,mm)-P.h1(i) );
        end
    end
    % store previous potentials:
    P.h2                    = O.h22(:,P.j,mm);
    % compute potentials at new cap:
    O.h22(:,P.j,mm)         = fnsyst( P, W );

    P.niter                 = P.niter+1;
    % Ribaltare il pezzo sopra con quello sotto (niter>10 va
    % all'inizio del while, per cui se si verifica condizione non
    % esegue il blocco di calcoli che stà adesso sopra).
    % Migliorare soluzione, o comunque salvare in M.nnc l'eventuale simulazione saltata
    if P.niter>10
        P.SS                = 1;
        if W.dt==W.dtmin
            P.LL            = 0;
        elseif W.dt>W.dtmin
            P.LL            = 1;
            W.dt            = W.dt/3;
            P.time(P.j)     = P.time(P.j)-2*W.dt;
            if W.dt<=W.dtmin
                P.time(P.j) = P.time(P.j) - W.dt + W.dtmin;
                W.dt        = W.dtmin;
            end
            % Ripristino flussi o potenziali nel caso di
            % annullamento dell'iterazione per riduzione del
            % W.dt.
            % Serve solo nel caso di P.LL=1, quando il ripristino
            % non può avvenire con le righe uguali poste alla fine
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
    W.dt                    = W.multmin * W.dt;
    if W.dt>W.dtmax
        W.dt                = W.dtmax;
    end
    P.LL                    = 0;

elseif and(P.niter>=4,P.niter<=10)
    W.dt                    = W.multmax * W.dt;
    if W.dt<W.dtmin
        W.dt                = W.dtmin;
    end
    P.LL                    = 0;
end