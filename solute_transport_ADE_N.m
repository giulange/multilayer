function P = solute_transport_ADE_N( P, W, S, B )
% P = solute_transport_ADE_N( P, W, S, B )
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

%% update intial concentrations

if P.j==1% first iteration
    % [W.dz,1:2,P.Nj]
    P.C1(:,1,P.j)       = S.CDE.Cin.NH;
    P.C1(:,2,P.j)       = S.CDE.Cin.NO;
elseif P.j>P.jstar% when?
    P.C1(:,:,P.j)       = P.C2(:,:,P.j-1); % --> assgnment before use? 
elseif P.j==P.jstar% when?
    P.C1(:,:,P.j)       = P.C1star(:,:,P.j);
end
P.C1star(:,:,P.j)       = P.C1(:,:,P.j);

% parte adsorbita
if W.ads==1
    % concentrazione in fase adsorbita (S Freundlich):
    P.S1(:,1,P.j)       = S.CDE.NX.Kf1(1)*P.C1(:,1,P.j).^S.CDE.NX.Kf2(1);
    P.S1(:,2,P.j)       = S.CDE.NX.Kf1(2)*P.C1(:,2,P.j).^S.CDE.NX.Kf2(2);
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
    idL                 = B.Ctop.dL/P.dz(1) +1; % <-- viene decimale, corretto??
    
    % Calcolo del fattore di riduzione del coefficiente di nitrificazione
    % (vedi Gusman et al, 1999. Analytical Modeling of Nitrogen Dynamics in
    % soil and groundwater. J.of Irrigation and Drainage Engineering.
    for k=1:idL
        if P.teta(k,P.j)<=P.tetafc(1)
            Fmw         = Fmw + P.teta(k,P.j)/P.tetafc(1); % / idL
        else
            Fmw         = Fmw + P.tetafc(1)/P.teta(k,P.j); % / idL
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
    % NH3 contribuisce per lo più nel momento della
    % somministrazione (per cui dopo non fai la sum):
    CNH3_UR             = CNH3_UR       + B.Ctop.Cstar.UR(i)*B.Ctop.KvUR*exp(-B.Ctop.KvUR*P.tm);
    CNH4_ORG_rp         = CNH4_ORG_rp   + B.Ctop.Cstar.ORG.rp(i)*B.Ctop.KmORG_rp*exp(-B.Ctop.KmORG_rp*tm);
    CNH4_ORG_sw         = CNH4_ORG_sw   + B.Ctop.Cstar.ORG.sw(i)*B.Ctop.KmORG_sw*exp(-B.Ctop.KmORG_sw*tm);

    % ?? DELETE ??
%     C_UR_RES(i)         = B.Ctop.Cstar.UR(i) - CNH4_UR_CUM(i) - CNH3_UR_CUM(i);
%     C_ORG_RES(i)        = B.Ctop.Cstar.ORG.rp(i) + B.Ctop.Cstar.ORG.sw(i) - CNH4_ORG_CUM(i);
end

% Non si moltiplica per W.dt perché dalla deirvata calcolata al for
% precedente si ottiene il CNH4 prodotto per unità di tempo che è quello
% che deve entrare nella equazione ADE.

% L'apporto della forma solida NH4 ed NO3 dura per l'intero periodo
% P.kk:P.kk+1. essendo in g/cm2/d, l'input è già nella forma richiesta
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
    fluxin              = P.fluxsurf(P.j);
    C1top               = P.Cinput(sl);
    lambdatop           = S.CDE.lambda(1);

    % --Non li usi--
%     Knitr_top           = S.CDE.Knitr(1);
%     Kimmb_top           = S.CDE.Kimmob(1);
%     Kdntr_top           = S.CDE.Kdenitr(i);
    % --Non li usi--

    while jj<=W.nlay
        for i=P.istar(jj)+1:P.istar(jj+1)
            % bot-cond for current layer
            % ...there is something not good in these statements...
            if i == P.istar(jj+1)
                P.C1(P.istar(jj+1)+1,sl,P.j)    = P.C1(P.istar(jj+1),sl,P.j);
                S.CDE.lambda(P.istar(jj+1)+1)   = S.CDE.lambda(P.istar(jj+1));
                S.CDE.Knitr(P.istar(jj+1)+1)    = S.CDE.Knitr(P.istar(jj+1));
                S.CDE.Kimmob(P.istar(jj+1)+1)   = S.CDE.Kimmob(P.istar(jj+1));
                S.CDE.Kdenitr(P.istar(jj+1)+1)  = S.CDE.Kdenitr(P.istar(jj+1));
            end

            % ?? DELETE ??
            % questa roba qui non la puoi scrivere qua dentro, semmai
            % andrebbe fuori di questa funzione!!
            P.dz(P.istar(jj+1))     = W.zint(jj+1)-P.z(P.istar(jj+1));
            P.dz(W.nz+1)            = P.dzbot;
            % ?? DELETE ??
            
            P.flux(W.nz+1,P.j)      = P.fluxbot(P.j);
            lambdaup                = (S.CDE.lambda(1)+lambdatop)/2;
            lambdalow               = (S.CDE.lambda(1)+S.CDE.lambda(2))/2;
            if i==P.isar(jj)+1
                fluxup              = ((P.flux(P.istar(jj)+1,P.j)+fluxin)/2);
                fluxlow             = ((P.flux(P.istar(jj)+1,P.j)+P.flux(P.istar(jj)+2,P.j))/2);
                C1up                = (P.C1(P.istar(jj)+1,sl,P.j)+C1top)/2;         
                C1low               = (P.C1(P.istar(jj)+1,sl,P.j)+P.C1(P.istar(jj)+2,sl,P.j))/2;
                if i>1
                    lambdaup        = (S.CDE.lambda(P.istar(jj)+1)+S.CDE.lambda(P.istar(jj)))/2;                
                    lambdalow       = (S.CDE.lambda(P.istar(jj)+1)+S.CDE.lambda(P.istar(jj)+2))/2;
                end
            else
                fluxup              = ((P.flux(i,P.j)+P.flux(i-1,P.j))/2);
                fluxlow             = ((P.flux(i,P.j)+P.flux(i+1,P.j))/2);
                C1up                = (P.C1(i,sl,P.j)+P.C1(i-1,sl,P.j))/2;                
                C1low               = (P.C1(i,sl,P.j)+P.C1(i+1,sl,P.j))/2;
                lambdaup            = (S.CDE.lambda(i) + S.CDE.lambda(i-1))/2;                
                lambdalow           = (S.CDE.lambda(i) + S.CDE.lambda(i+1))/2;
            end

            % The following is the discretised Eq. for Cphys (tracer):
            % ...
            C2one                   = (fluxup*C1up - fluxlow*C1low)/P.dz(i);
            if i==P.istar(jj)+1
                C2two               = lambdaup*abs(fluxup)*(C1top-P.C1(i,sl,P.j)) / P.dztop;
            else                
                C2two               = lambdaup*abs(fluxup)*(P.C1(i-1,sl,P.j)-P.C1(i,sl,P.j)) / P.dz(i);
            end
            C2three                 = lambdalow*abs(fluxlow)*(P.C1(i,sl,P.j)-P.C1(i+1,sl,P.j))/P.dz(i);
            C2_phys                 = W.dt*(-C2one + ((C2two-C2three)/P.dz(i)));
            % [g cm-3 H2O]
            C2phys                  = ( C2_phys + P.C1(i,sl,P.j)*P.teta(i,P.j) ) / P.teta(i,P.j);
            
            if P.teta(i,P.j)<P.tetafc(i)
                teta_ratio          = P.teta(i,P.j)/P.tetafc(i);
            else
                teta_ratio          = P.tetafc(i)/P.teta(i,P.j);
            end
            C2_ntf_lq               = S.CDE.Knitr(i)*1.07^(B.Ctop.Tstar(P.kk)-S.CDE.Topt)*teta_ratio*P.C1(i,1,P.j)*P.teta(i,P.j);
            C2_ntf_sd               = S.CDE.Knitr(i)*1.07^(B.Ctop.Tstar(P.kk)-S.CDE.Topt)*teta_ratio*P.dap(i)*P.S1(i,1,P.j);
            % attingimento selettivo:
            C2_sink                 = S.CDE.NX.Kr(sl)*P.sink(i,P.j)*P.C1(i,sl,P.j);

            if sl==1
                % SsSk: Source-Sink for solutes [g cm-3 suolo day-1]
                SsSk                = -C2_ntf_lq -C2_ntf_sd -C2_sink + CNH4_pn(i);
            elseif sl==2
                % Sm(z,t), Eq. 16
                C2_immb             = S.CDE.Kimmb(i)*1.05^(B.Ctop.Tstar(P.kk)-S.CDE.Topt)*teta_ratio*P.C1(i,2,P.j)*P.teta(i,P.j);
                
                % teta_tsh: the threshold water content for de-nitrification:
                teta_tsh            = 0.627*P.tetafc(i)-0.0267*(P.tetas(i)-P.teta(i,P.j))/P.tetas(i)*P.tetafc(i);
                if P.teta(i,P.j)>teta_tsh
                    % Sd(z,t), Eq. 17
                    C2_dntf         = S.CDE.Kdntr(i)*1.07^(B.Ctop.Tstar(P.kk)-S.CDE.Topt)*(P.teta(i,P.j)-teta_tsh)/(P.tetafc(i)-teta_tsh)*P.C1(i,2,P.j)*P.teta(i,P.j);
                else
                    C2_dntf         = 0;
                end
                % SsSk: Source-Sink for solutes [g cm-3 suolo day-1]
                SsSk                = +C2_ntf_lq +C2_ntf_sd -C2_immb -C2_dntf -C2_sink +CNO_pn(i);
            end
            
            % è giusto aggiungere P.C1?? Risponde Antonio...
            C2_chem                 = ( SsSk*W.dt + P.C1(i,sl,P.j)*P.teta(i,P.j) ) /P.teta(i,P.j);
            C2_tot                  = C2phys + C2_chem;

            if S.CDE.NX.Kf2(sl)==1
                P.C2(i,sl,P.j)      = ( C2phys*P.teta(i,P.j) + P.dap(i)*S.CDE.NX.Kf1(sl)*P.C1(i,sl,P.j) + SsSk*W.dt ) ...
                                              / ( fnteta( W.hfc, P, i ) + P.dap(i)*S.CDE.NX.Kf1(sl) );
            else    
%% Bisection method
% Bisection method (vedi libro con esercizi numerici in Matlab) for
% adsorbing solutes with Freundlich f_bis_sol è la funzione di cui si vuole
% trovare lo zero.
% E' bene fissare un limite inferiore per P.C2(i,P.j)=10^k per eliminare
% eventuali valori negativi ed evitare problemi di convergenza del metodo.
                if C2_tot < 10^-9
                    C2phys          = 0;
                    P.C2(i,sl,P.j)  = C2phys;
                else
                    P.C2(i,sl,P.j)  = mln_bisection( P.teta(i,P.j),P.dap(i),P.S1(i,sl,P.j), ...
                                        S.CDE.NX.Kf1(sl),S.CDE.NX.Kf2(sl),C2phys,SsSk,W.dt);
%                     f_bis_sol       = @(yps) P.teta(i,P.j)*yps-P.teta(i,P.j)*C2phys+P.dap(i)*S.CDE.NX.Kf1(sl)*yps^S.CDE.NX.Kf2(sl)-P.dap(i)*P.S1(i,sl,P.j)-SsSk*W.dt;
%                     a_bs            = -100;
%                     b_bs            = 100;
%                     espon           = -30;
%                     eps_bs          = 10^espon;
%                     fa              = feval(f_bis_sol,a_bs);
%                     fb              = feval(f_bis_sol,b_bs);
%                     if fa*fb>0
%                         error('Intervallo iniziale non accettabile')
%                     end
%                     Nit_sol         = (log(b_bs-a_bs)-log(eps_bs))/log(2);
%                     for u=3:Nit_sol+2
%                         ics         = (a_bs+b_bs)/2;
%                         fics        = feval(f_bis_sol,ics);
%                         if fa*fics<=0
%                             b_bs    = ics;
%                         else
%                             a_bs    = ics;
%                             fa      = fics;
%                         end
%                     end
%                     ics_fin         = (a_bs+b_bs)/2;
%                     P.C2(i,sl,P.j)  = ics_fin;
                end
            end
        end
        %-----------------------------------------------------------------
        % operazioni prima di passare al prossimo strato:
        P.C2(P.istar(jj+1)+1,sl,P.j)        = P.C2(P.istar(jj+1),sl,P.j);

        % USELESS --> ??DELETE??
        C1top                               = P.C2(P.istar(jj+1)+1,sl,P.j);
        fluxin                              = P.flux(P.istar(jj+1)+1);
        % USELESS --> ??DELETE??
        
        jj                                  = jj+1;
        %-----------------------------------------------------------------
    end
end
% -----------------------------------------------------------------------
% soluzione equazione convezione dispersione :: END
% -----------------------------------------------------------------------
return