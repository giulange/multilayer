%% multilayer_wtus_thomas
% Thomas algorithm, using previous approach implemented by Antonio Coppola.
%% lettura h iniziali o aggiornamento h -- check at the end!
% lettura h iniziali e aggiornamento h
if P.j==1
    P.h1                    = P.hin;
else %if P.j>P.jstar  % if simulation runs normally
    P.h1                    = O.h22(:,P.j-1); % --> not yet initialized!!
end
%% teta, cond, cap, sink :: all nodes
% calc. param. ritenz. e conducib. ai vari nodi nei diversi strati
for i=1:P.nz
    % Ritenzione:
    P.teta(i)               = fnteta( P.h1, P.sh, i );
    % Conducibilità:
    if P.sh.ifc(i)==1 || P.sh.ifc(i)==3
        P.kond(i)           = multilayer_conductivity_node( P.teta,P.sh, i );
    else
        P.kond(i)           = multilayer_conductivity_node( P.h1,  P.sh, i );
    end
    % Capacità:
    P.cap(i)                = multilayer_capacity(  P.h1,  P.sh, i );
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
% per infiltrazione:
P.fluxsurf_max(P.j)         = -P.km_max(P.j)*((W.hsurfmax-P.h1(1))/(P.nodes.dz(1)/2)+1);% Darcy
% ponding --> W.hsurfmax = [0,5] cm, pressione che si sviluppa nel nodo
% fittizio.

switch W.itbc% 0:flux; 1:potential;
    case 0
        if and( W.qsurf==0, W.iveg==1 )
            W.hsurf         = 13.3 * 10^5 * log(W.vpr);
            P.teta_hsurf    = fnteta( W.hsurf, P.sh, 1 );

            if or(P.sh.ifc==1,P.sh.ifc==3)
%                 P.km(P.j)   = (P.kond(1) + multilayer_conductivity_node(P.teta_hsurf,P.sh,1))/2;
                P.km(P.j)   = multilayer_conductivity_internode( ...
                                multilayer_conductivity_node(P.teta_hsurf,P.sh,1), ...
                                P.kond(1), W.Kmeth );
            else
%                 P.km(P.j)   = (P.kond(1) + multilayer_conductivity_node(W.hsurf,P.sh,1))/2;
                P.km(P.j)   = multilayer_conductivity_internode( ...
                                multilayer_conductivity_node(W.hsurf,P.sh,1), ...
                                P.kond(1), W.Kmeth );
            end
            P.Emax          = -P.km(P.j) * ((W.hsurf-P.h1(1))/(P.nodes.dz(1)/2)+1);

            if P.Ep<P.Emax
                W.qsurf     = P.Ep;
                W.hsurf     = (3*P.h1(1)-P.h1(2))/2;
                O.fluxsurf(1,P.j,mm) = W.qsurf; % porta fuori (ma prima controlla che sia qsurf e non hsurf!!
            else
                W.qsurf     = P.Emax;
                O.fluxsurf(1,P.j,mm) = W.qsurf; % porta fuori
            end
            
        elseif W.qsurf<0
            W.hsurf         = (3*P.h1(1)-P.h1(2))/2;
            O.fluxsurf(1,P.j,mm) = W.qsurf;
        elseif and(W.qsurf==0,W.iveg==0)
            W.hsurf         = (3*P.h1(1)-P.h1(2))/2;
            O.fluxsurf(1,P.j,mm) = W.qsurf;
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
            P.rnf           = 1;
            W.itbc          = 1;
        end

    case 1
        if P.rnf==1
            W.hsurf         = W.hsurfmax;
            P.km_max(P.j)   = (P.kond(1)+P.sh.k0(1))/2;
            P.fluxsurf_max(P.j)     = -P.km_max(P.j) * ((W.hsurfmax-P.h1(1))/(P.nodes.dz(1)/2)+1);
            O.fluxsurf(1,P.j,mm)    = P.fluxsurf_max(P.j);
        elseif P.rnf==0
            P.teta_hsurf    = fnteta( W.hsurf, P.sh, 1 );
            if or(P.sh.ifc==1,P.sh.ifc==3)
%                 P.km(P.j)   = ( P.kond(1) + multilayer_conductivity_node(P.teta_hsurf,P.sh,1) )/2;
                P.km(P.j)   = multilayer_conductivity_internode( ...
                                multilayer_conductivity_node(P.teta_hsurf,P.sh,1), ...
                                P.kond(1), W.Kmeth );
            else
%                 P.km(P.j)   = ( P.kond(1) + multilayer_conductivity_node(W.hsurf,P.sh,1) )/2;
                P.km(P.j)   = multilayer_conductivity_internode( ...
                                multilayer_conductivity_node(W.hsurf,P.sh,1), ...
                                P.kond(1), W.Kmeth );
            end
            O.fluxsurf(1,P.j,mm)    = -P.km(P.j) * ((W.hsurf-P.h1(1))/(P.nodes.dz(1)/2)+1);
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
    O.fluxbot(1,P.j,mm)     = W.qbot;
    W.hbot(P.j)             = (P.h1(P.nz)-P.h1(P.nz-1))*...
                               P.nodes.dz(P.nz+1)/P.nodes.dz(P.nz)+P.h1(P.nz);
%             case {1,2}
elseif or(W.ibbc==1,W.ibbc==2)
    P.teta_hbot             = fnteta( W.hbot, P.sh, P.nz);
    if or(P.sh.ifc==1,P.sh.ifc==3)
%         P.kp(P.nz)          = ( P.kond(P.nz) + multilayer_conductivity_node(P.teta_hbot,P.sh,P.nz) )/2;
        P.kp(P.nz)          = multilayer_conductivity_internode( P.kond(P.nz), ...
                                multilayer_conductivity_node(P.teta_hbot,P.sh,P.nz), W.Kmeth );
    else
%         P.kp(P.nz)          = ( P.kond(P.nz) + multilayer_conductivity_node(W.hbot,P.sh,P.nz) )/2;
        P.kp(P.nz)          = multilayer_conductivity_internode( P.kond(P.nz), ...
                                multilayer_conductivity_node(W.hbot,P.sh,P.nz), W.Kmeth );
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
%% calcolo flussi all'interfaccia -- check with Antonio (+2?)
% Definiamo l'obiettivo, poi valutiamo come implementare.
% P.k = 2:W.nlay;

% why +2 and not +1??
P.flux(P.nodes.cumsum(2:W.nlay)) = 2/3 * P.flux(P.nodes.cumsum(2:W.nlay)-1) +...
                                   1/3 * P.flux(P.nodes.cumsum(2:W.nlay)+2);
%% flussi cumulati in superficie ==> not used
%         if P.j==1
%             fluxsurfcum(P.j)        = O.fluxsurf(1,P.j,mm)*P.dt;
%         else
%             fluxsurfcum(P.j)        = O.fluxsurf(1,P.j,mm)*P.dt + ...
%                                       fluxsurfcum(P.j-1);
%         end
%% flussi cumulati ai nodi intermedi ==> not used -- error on P.j
%         for i=1:P.nz
%             if P.j==1
%                 fluxcum(i,1)        = P.flux(i,1)*P.dt;
%             else
%                 fluxcum(i,P.j)      = P.flux(i,P.j)*P.dt + ...
%                                       fluxcum(i,P.j-1);
%             end
%         end
%% soluzione del sistema tridiagonale (see WARNING) -- obiettivo?
% Obiettivo? h al j+1
P.h2                        = fnsyst( P, W );% al j+1
%% calcolo capacita' metodo implicito (vedi Karvonen) -- capacità 1/2
for i=1:P.nz
    if abs(P.h2(i)-P.h1(i))>W.tolle1
        P.cap(i)            = ( fnteta(P.h2,P.sh,i) - ...
                            fnteta(P.h1,P.sh,i) ) / ( P.h2(i)-P.h1(i) );
    end
end
% fnsyst con il nuovo P.cap:
O.h22(:,P.j,mm)             = fnsyst( P, W );
%         figure(13),hold on, plot(squeeze(O.h22(:,P.j,mm))), hold off;
%% calcolo capacita' al tempo medio -- obiettivo?
% W.tolle1=tolleranza per la convergenza della capacità
% P.niter                     = 0; % Nj per convergenza capacità.
% while max(abs(O.h22(:,P.j,mm)-P.h2))>W.tolle1
for p = 0:W.maxit
    
    if max(abs(O.h22(:,P.j,mm)-P.h2))<W.tolle1
        nr_breaked  = true;
        break
    end

    for i=1:P.nz            
        if abs(O.h22(i,P.j,mm)-P.h2(i))>W.tolle1        
            P.cap(i)        = ( fnteta(O.h22(:,P.j,mm),P.sh,i) - ...
                                fnteta(P.h1,P.sh,i) ) / ( O.h22(i,P.j,mm)-P.h1(i) );
        end
    end
    % store previous potentials:
    P.h2                    = O.h22(:,P.j,mm);
    % compute potentials at new cap:
    O.h22(:,P.j,mm)         = fnsyst( P, W );

%     niter                 = niter+1;
%     if niter > W.maxit, break; end
end

if nr_breaked
    P.iter(:,P.j) = [p;nr_breaked;fl_noconv;n_noconv];
    % Copied by HYDRUS(?):
    if p<=3
        P.dt                = W.multmin * P.dt;
        if P.dt>W.dtmax
            P.dt            = W.dtmax;
        end
    elseif and(p>=W.maxit-3,p<=W.maxit)
        P.dt                = W.multmax * P.dt;
        if P.dt<W.dtmin
            P.dt            = W.dtmin;
        end
    end
else
    % Migliorare soluzione, o comunque salvare in M.nnc l'eventuale
    % simulazione saltata
    if P.dt==W.dtmin
        % *CONTINUE WITHOUT CONVERGENCE:
        fl_noconv = true;% is non-convergent?
        % ...what to do next? --> we should do something
    elseif P.dt>W.dtmin
        n_noconv = n_noconv +1;
        P.iter(:,P.j) = [p;nr_breaked;fl_noconv;n_noconv];
        P.dt            = P.dt/3;
        P.time(P.j)     = P.time(P.j)-2*P.dt;
        if P.dt<=W.dtmin
            P.time(P.j) = P.time(P.j) - P.dt + W.dtmin;
            P.dt        = W.dtmin;
        end
        
        % If we go back to the previous integer, we must update:
        P.tidx              = floor(P.time(P.j))+1;

% --- useless        
        % Ripristino flussi o potenziali nel caso di
        % annullamento dell'iterazione per riduzione del
        % P.dt.
        % Serve solo nel caso di P.LL=1, quando il ripristino
        % non può avvenire con le righe uguali poste alla fine
        % del blocco per la simulazione del trasporto di
        % soluti.
%         if W.itbc==0
%             W.qsurf     = B.top.hqstar(P.tidx);
%         elseif and(W.itbc==1,P.rnf==0)
%             W.hsurf     = B.top.hqstar(P.tidx);
%         elseif and(W.itbc==1,P.rnf==1)
%             W.qsurf     = B.top.hqstar(P.tidx);
%         end
        % Ripristino tempo di lettura delle EC.data nel caso in
        % cui l'annullamento della simulazione avvenga a
        % cavallo del cambio di EC.t
        if W.iosm==1 && V.ifs>3 && P.IEC==1
            P.ktec      = P.ktec-1;
        end
% --- useless
        
    end
end
