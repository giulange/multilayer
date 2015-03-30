function P = solute_transport_ADE_N( P, W, S, B, O, mm )
% P = solute_transport_ADE_N( P, W, S, B, O, mm )
% 
% DESCRIPTION
%   This function accounts for two forms of N, namely NO3 and NH4.
%   This module accounts for the following processes [see the paper by
%   Zue_ENVIRONMENTAL POLLUTION_2007]:
%       UREA
%               --> IDROLISI            (pedice 'h')    --> NH4
%               --> VOLATILIZZAZIONE    (pedice 'v')    --> NH3
%       MATERIA ORGANICA
%               --> MINERALIZZAZIONE    (pedice 'm')    --> NH4
% 
%   For nitrification, immobilization and denitrification see She et al.,
%   2005. An Inverse Method to Estimate the Source-Sink Term in the Nitrate
%   Transport Equation. SSSAJ Soil Sci. Soc. Am. J. 71:26-34 
% 
%   Remember that sl accounts for the type of solute:
%       sl=1    --> NH_FR
%       sl=2    --> NO_FR
%   while P.kk accounts for Ctopboundary input (see B.Ctop).
% 
% 
% INPUTS
%   P:      The structure array with all parameters generated during the
%           current simulation.
%   W:      ??
%   S:      ??
%   B:      ??
% 
% 
% OUTPUTS
%   P:      Newly created variables are stored in P and returned to the
%           main function of the program. (which variables?)

%% init
% 
jj                      = 1;
Fmw                     = 0;
CNH4_UR                 = 0;
CNH3_UR                 = 0; %--> la calcoli ma non la usi!!!!
CNH4_ORG_rp             = 0;
CNH4_ORG_sw             = 0;
CNH4_pn                 = NaN(W.nz,1);
CNO_pn                  = NaN(W.nz,1);
C1                      = NaN(W.nz,2);
S1                      = NaN(W.nz,2);
%% update intial concentrations

if P.j==1% first iteration
    % [W.dz,1:2,P.Nj]
    C1(:,1)             = S.CDE.Cin.NH;
    C1(:,2)             = S.CDE.Cin.NO;
elseif P.j>P.jstar% se avanza...
    C1                  = O.C2(:,P.j-1,mm,:); % --> assgnment before use? 
elsef P.j==P.jstar% se non puo' avanzare...
    C1                  = O.C2(:,P.j-1,mm,:);
end

% parte adsorbita
if W.ads==1
    % concentrazione in fase adsorbita (S Freundlich):
    S1(:,1)       = S.CDE.NX.Kf1(1)*C1(:,1).^S.CDE.NX.Kf2(1);
    S1(:,2)       = S.CDE.NX.Kf1(2)*C1(:,2).^S.CDE.NX.Kf2(2);
end
%% ?? define cell / delete cell ??                

% ?? aggiustare su (sl) ??
if W.iCtopvar==0
    if and(P.time(P.j)>=S.CDE.tCinput,P.CC==0)
        P.Cinput        = S.CDE.Cinput;
        if P.time(P.j)>=S.CDE.tCinput_end
            P.CC        = 1;
        end
    else
        P.Cinput        = 0;
    end
end
%% IDROLISI E VOLATILIZZAZIONE UREA + MINERALIZZAZIONE SOM

% -----------------------------------------------------------------------
% mineralizzazione sostanza organica :: START
% -----------------------------------------------------------------------
for i=1:P.kk
    tm                  = P.time(P.j) - B.top.thqstar(i);

%     P.CNH4_UR_CUM(i) = B.Ctop.Cstar.UR(i)*(1-exp(-B.Ctop.KhUR*P.tm));
%     P.CNH3_UR_CUM(i) = B.Ctop.Cstar.UR(i)*(1-exp(-B.Ctop.KvUR*P.tm));      

    % Calcolo numero di nodi nello strato di interramento
    % del concime azotato:
    idL                 = B.Ctop.dL/P.nodes.dz(1) +1; % <-- viene decimale, corretto??
    
    % Calcolo del fattore di riduzione del coefficiente di nitrificazione
    % (vedi Gusman et al, 1999. Analytical Modeling of Nitrogen Dynamics in
    % soil and groundwater. J.of Irrigation and Drainage Engineering.
    for k=1:idL
        if P.teta(k,P.j)<=P.sh.tetafc(1)
            Fmw         = Fmw + P.teta(k,P.j)/P.sh.tetafc(1); % / idL
        else
            Fmw         = Fmw + P.sh.tetafc(1)/P.teta(k,P.j); % / idL
        end
    end
    % Il fattore di riduzione non viene calcolato nodo per nodo ma viene
    % mediato sull'intero strato di interramento.
    Fmw                 = Fmw / idL;

    B.Ctop.KmORG_rp     = 5.6*10^12 * exp(-9800/(B.Ctop.Tstar(P.kk)+273))*Fmw;
    B.Ctop.KmORG_sw     = 4.0*10^9  * exp(-8400/(B.Ctop.Tstar(P.kk)+273))*Fmw;

%     if and(P.time(P.j)>40,P.time(P.j)<70)
%         B.Ctop.KmORG_sw = 4.0*10^9.5 * ...
%             exp(-8400/(B.Ctop.Tstar(P.kk)+273))*P.Fmw;
%     end
%     if P.time(P.j)>95
%         B.Ctop.KmORG_sw = 4.0*10^3.0 * ...
%             exp(-8400/(B.Ctop.Tstar(P.kk)+273))*P.Fmw;
%     end
% 
%     P.CNH4_ORG_rp_CUM(i) = B.Ctop.Cstar.ORG.rp(i) * ...
%                             (1-exp(-B.Ctop.KmORG_rp*P.tm));
%     P.CNH4_ORG_sw_CUM(i) = B.Ctop.Cstar.ORG.sw(i) * ...
%                             (1-exp(-B.Ctop.KmORG_sw*P.tm));
%     P.CNH4_ORG_CUM(i) = P.CNH4_ORG_rp_CUM(i)  + ...
%                             P.CNH4_ORG_sw_CUM(i);

    % qui applica il decadimento (asintotico, quanto viene prodotto):
    CNH4_UR             = CNH4_UR       + B.Ctop.Cstar.UR(i)*B.Ctop.KhUR*exp(-B.Ctop.KhUR*tm);
    % NH3 contribuisce per lo pi� nel momento della
    % somministrazione (per cui dopo non fai la sum):
    CNH3_UR             = CNH3_UR       + B.Ctop.Cstar.UR(i)*B.Ctop.KvUR*exp(-B.Ctop.KvUR*tm);
    CNH4_ORG_rp         = CNH4_ORG_rp   + B.Ctop.Cstar.ORG.rp(i)*B.Ctop.KmORG_rp*exp(-B.Ctop.KmORG_rp*tm);
    CNH4_ORG_sw         = CNH4_ORG_sw   + B.Ctop.Cstar.ORG.sw(i)*B.Ctop.KmORG_sw*exp(-B.Ctop.KmORG_sw*tm);

    % ?? DELETE ??
%     C_UR_RES(i)         = B.Ctop.Cstar.UR(i) - CNH4_UR_CUM(i) - CNH3_UR_CUM(i);
%     C_ORG_RES(i)        = B.Ctop.Cstar.ORG.rp(i) + B.Ctop.Cstar.ORG.sw(i) - CNH4_ORG_CUM(i);
end

% Non si moltiplica per W.dt perch� dalla deirvata calcolata al for
% precedente si ottiene il CNH4 prodotto per unit� di tempo che � quello
% che deve entrare nella equazione ADE.

% L'apporto della forma solida NH4 ed NO3 dura per l'intero periodo
% P.kk:P.kk+1. essendo in g/cm2/d, l'input � gi� nella forma richiesta
% dall'ADE.
% Questo viene poi moltiplicato per W.dt nella stessa equazione.
CNH4_SD                 = B.Ctop.Cstar.NH.SD(P.kk);
CNO_SD                  = B.Ctop.Cstar.NO.SD(P.kk);
% Costruzione dei pool:
POOL_NH4_SD             = CNH4_UR + CNH4_ORG_rp + CNH4_ORG_sw + CNH4_SD;
POOL_NO_SD              = CNO_SD; % [g cm-2]
% Una volta costituito, il POOL si assume distribuito per l'intero spessore
% di suolo B.Ctop.dL ==> lo si divide per B.Ctop.dL e si ottiene una
% concnetrazione in g/cm3 di suolo.
for i=1:W.nz
    if i<=idL
        % se sono nello spessore di interramento:
        CNH4_pn(i)      = POOL_NH4_SD / B.Ctop.dL;
        CNO_pn(i)       = POOL_NO_SD  / B.Ctop.dL;
    else
        % se sono sotto lo spessore di interramento:
        CNH4_pn(i)      = 0;
        CNO_pn(i)       = 0;
    end
end
% -----------------------------------------------------------------------
% mineralizzazione sostanza organica :: END
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
% soluzione equazione convezione dispersione :: START
% -----------------------------------------------------------------------
for sl=1:2
    fluxin              = O.fluxsurf(1,P.j,mm);
    C1top               = P.Cinput(sl);
    lambdatop           = S.CDE.lambda(1);

    % --Non li usi--
%     Knitr_top           = S.CDE.Knitr(1);
%     Kimmb_top           = S.CDE.Kimmob(1);
%     Kdntr_top           = S.CDE.Kdenitr(i);
    % --Non li usi--

    while jj<=W.nlay                                    % |--> substitute these two with "for inode = 1:W.nz, ..., end"
        for i=P.nodes.cumsum(jj)+1:P.nodes.cumsum(jj+1) %/
            % bot-cond for current layer
            % ...there is something not good in these statements...
            if i == P.nodes.cumsum(jj+1)
                C1(P.nodes.cumsum(jj+1)+1,sl)          = C1(P.nodes.cumsum(jj+1),sl);
                S.CDE.lambda(P.nodes.cumsum(jj+1)+1)   = S.CDE.lambda(P.nodes.cumsum(jj+1));
                S.CDE.Knitr(P.nodes.cumsum(jj+1)+1)    = S.CDE.Knitr(P.nodes.cumsum(jj+1));
                S.CDE.Kimmob(P.nodes.cumsum(jj+1)+1)   = S.CDE.Kimmob(P.nodes.cumsum(jj+1));
                S.CDE.Kdenitr(P.nodes.cumsum(jj+1)+1)  = S.CDE.Kdenitr(P.nodes.cumsum(jj+1));
            end

            P.flux(W.nz+1,P.j)      = O.fluxbot(1,P.j,mm);
            lambdaup                = (S.CDE.lambda(1)+lambdatop)/2;
            lambdalow               = (S.CDE.lambda(1)+S.CDE.lambda(2))/2;
            if i==P.nodes.cumsum(jj)+1
                fluxup              = ((P.flux(P.nodes.cumsum(jj)+1,P.j)+fluxin)/2);
                fluxlow             = ((P.flux(P.nodes.cumsum(jj)+1,P.j)+P.flux(P.nodes.cumsum(jj)+2,P.j))/2);
                C1up                = (C1(P.nodes.cumsum(jj)+1,sl)+C1top)/2;         
                C1low               = (C1(P.nodes.cumsum(jj)+1,sl)+C1(P.nodes.cumsum(jj)+2,sl))/2;
                if i>1
                    lambdaup        = (S.CDE.lambda(P.nodes.cumsum(jj)+1)+S.CDE.lambda(P.nodes.cumsum(jj)))/2;                
                    lambdalow       = (S.CDE.lambda(P.nodes.cumsum(jj)+1)+S.CDE.lambda(P.nodes.cumsum(jj)+2))/2;
                end
            else
                fluxup              = ((P.flux(i,P.j)+P.flux(i-1,P.j))/2);
                fluxlow             = ((P.flux(i,P.j)+P.flux(i+1,P.j))/2);
                C1up                = (C1(i,sl)+C1(i-1,sl))/2;                
                C1low               = (C1(i,sl)+C1(i+1,sl))/2;
                lambdaup            = (S.CDE.lambda(i) + S.CDE.lambda(i-1))/2;                
                lambdalow           = (S.CDE.lambda(i) + S.CDE.lambda(i+1))/2;
            end

            % The following is the discretised Eq. for Cphys (tracer):
            % ...
            C2one                   = (fluxup*C1up - fluxlow*C1low)/P.nodes.dz(i);
            if i==P.nodes.cumsum(jj)+1
                C2two               = lambdaup*abs(fluxup)*(C1top-C1(i,sl)) / P.nodes.dz(1);
            else                
                C2two               = lambdaup*abs(fluxup)*(C1(i-1,sl)-C1(i,sl)) / P.nodes.dz(i);
            end
            C2three                 = lambdalow*abs(fluxlow)*(C1(i,sl)-C1(i+1,sl))/P.nodes.dz(i);
            C2_phys                 = W.dt*(-C2one + ((C2two-C2three)/P.nodes.dz(i)));
            % [g cm-3 H2O]
            C2phys                  = ( C2_phys + C1(i,sl)*P.teta(i,P.j) ) / P.teta(i,P.j);
            
            if P.teta(i,P.j)<P.sh.tetafc(i)
                teta_ratio          = P.teta(i,P.j)/P.sh.tetafc(i);
            else
                teta_ratio          = P.sh.tetafc(i)/P.teta(i,P.j);
            end
            C2_ntf_lq               = S.CDE.Knitr(i)*1.07^(B.Ctop.Tstar(P.kk)-S.CDE.Topt)*teta_ratio*C1(i,1)*P.teta(i,P.j);
            C2_ntf_sd               = S.CDE.Knitr(i)*1.07^(B.Ctop.Tstar(P.kk)-S.CDE.Topt)*teta_ratio*P.sh.dap(i)*S1(i,1);
            % attingimento selettivo:
            C2_sink                 = S.CDE.NX.Kr(sl)*P.sink(i,P.j)*C1(i,sl);

            if sl==1
                % SsSk: Source-Sink for solutes [g cm-3 suolo day-1]
                SsSk                = -C2_ntf_lq -C2_ntf_sd -C2_sink + CNH4_pn(i);
            elseif sl==2
                % Sm(z,t), Eq. 16
                C2_immb             = S.CDE.Kimmb(i)*1.05^(B.Ctop.Tstar(P.kk)-S.CDE.Topt)*teta_ratio*C1(i,2)*P.teta(i,P.j);
                
                % teta_tsh: the threshold water content for de-nitrification:
                teta_tsh            = 0.627*P.sh.tetafc(i)-0.0267*(P.sh.tetas(i)-P.teta(i,P.j))/P.sh.tetas(i)*P.sh.tetafc(i);
                if P.teta(i,P.j)>teta_tsh
                    % Sd(z,t), Eq. 17
                    C2_dntf         = S.CDE.Kdntr(i)*1.07^(B.Ctop.Tstar(P.kk)-S.CDE.Topt)*(P.teta(i,P.j)-teta_tsh)/(P.sh.tetafc(i)-teta_tsh)*C1(i,2)*P.teta(i,P.j);
                else
                    C2_dntf         = 0;
                end
                % SsSk: Source-Sink for solutes [g cm-3 suolo day-1]
                SsSk                = +C2_ntf_lq +C2_ntf_sd -C2_immb -C2_dntf -C2_sink +CNO_pn(i);
            end
            
            % � giusto aggiungere P.C1?? Risponde Antonio...
            C2_chem                 = ( SsSk*W.dt + C1(i,sl)*P.teta(i,P.j) ) /P.teta(i,P.j);
            C2_tot                  = C2phys + C2_chem;

            if S.CDE.NX.Kf2(sl)==1
                O.C2(i,P.j,mm,sl) = ( C2phys*P.teta(i,P.j) + P.sh.dap(i)*S.CDE.NX.Kf1(sl)*C1(i,sl) + SsSk*W.dt ) ...
                                              / ( fnteta( W.hfc, P.sh, i ) + P.sh.dap(i)*S.CDE.NX.Kf1(sl) );
            else    
%% Bisection method
% Bisection method (vedi libro con esercizi numerici in Matlab) for
% adsorbing solutes with Freundlich f_bis_sol � la funzione di cui si vuole
% trovare lo zero.
% E' bene fissare un limite inferiore per O.C2(i,P.j,mm,sl)=10^k per eliminare
% eventuali valori negativi ed evitare problemi di convergenza del metodo.
                if C2_tot < 10^-9
                    C2phys          = 0;
                    O.C2(i,P.j,mm,sl) = C2phys;
                else
                f_bis_sol       = @(yps) ...
                    P.teta(i,P.j)*yps        - ...
                    P.teta(i,P.j)*C2phys     + ...
                    P.sh.dap(i)*S.CDE.NX.Kf1(sl)*yps^S.CDE.NX.Kf2(sl) - ...
                    P.sh.dap(i)*S1(i,sl)          - ...
                    SsSk*W.dt;

                    O.C2(i,P.j,mm,sl)  = mln_bisection( f_bis_sol, -100, 100, -30 );
                    %mln_bisection( P.teta(i,P.j),P.sh.dap(i),S1(i,sl), ...
                                        %S.CDE.NX.Kf1(sl),S.CDE.NX.Kf2(sl),C2phys,SsSk,W.dt );
                end
            end
        end
        %-----------------------------------------------------------------
        % operazioni prima di passare al prossimo strato:
        O.C2(P.nodes.cumsum(jj+1)+1,P.j,mm,sl)     = O.C2(P.nodes.cumsum(jj+1),P.j,mm,sl);

        % USELESS --> ??DELETE??
        C1top                               = O.C2(P.nodes.cumsum(jj+1)+1,P.j,mm,sl);
        fluxin                              = P.flux(P.nodes.cumsum(jj+1)+1);
        % USELESS --> ??DELETE??
        
        jj                                  = jj+1;
        %-----------------------------------------------------------------
    end
end
% -----------------------------------------------------------------------
% soluzione equazione convezione dispersione :: END
% -----------------------------------------------------------------------
return