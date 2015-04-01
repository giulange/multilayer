function [C2out,P] = solute_transport_ADE_N( P, W, S, B, C1, fluxsurf, fluxbot )
% [C2out,P] = solute_transport_ADE_N( P, W, S, B, C1, fluxin, fluxout )
% 
% OLD CALL:
%   [O,P] = solute_transport_ADE_N( P, W, S, B, O, mm )
% 
% NOTES
%   I should avoid passing a lot of variables. In particular I should use a
%   more compact O in terms of used dimensions; i.e. give in input
%   O.C2(:,[P.j-1:P.j],mm,:) as squeeze( O.C2(:,2,1,:) ) --> O.C2(:,2,:).
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
%   while P.kk accounts for Ctopboundary input (see B).
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
Fmw                     = 0;
CNH4_UR                 = 0;
CNH3_UR                 = 0; %--> la calcoli ma non la usi!!!!
CNH4_ORG_rp             = 0;
CNH4_ORG_sw             = 0;
CNH4_pn                 = NaN(P.nz,1);
CNO_pn                  = NaN(P.nz,1);
S1                      = NaN(P.nz,2);
C2out                   = NaN(P.nz,2);
%% update intial concentrations
% if P.j==1% first iteration
%     % [W.dz,1:2,P.Nj]
%     C1(:,1)             = S.Cin.NH;
%     C1(:,2)             = S.Cin.NO;
% elseif P.j>P.jstar% se avanza...
%     % OLD % P.C1(:,:,P.j) = P.C2(:,:,P.j-1);
%     C1                  = C2in;
% elseif P.j==P.jstar% se non puo' avanzare...
%     % OLD % P.C1(:,:,P.j) = P.C1star(:,:,P.j);
%     C1                  = C2in;
% end
%% parte adsorbita
if W.ads==1
    % concentrazione in fase adsorbita (S Freundlich):
    S1(:,1)       = S.NX.Kf1(1)*C1(:,1).^S.NX.Kf2(1);
    S1(:,2)       = S.NX.Kf1(2)*C1(:,2).^S.NX.Kf2(2);
end
%% ?? define cell / delete cell ??
% ?? aggiustare su (sl) ??
if W.iCtopvar==0
    if and(P.time(P.j)>=S.tCinput,P.CC==0)
        P.Cinput        = S.Cinput;
        if P.time(P.j)>=S.tCinput_end
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

% Calcolo numero di nodi nello strato di interramento
% del concime azotato:
%     idL                 = B.dL/P.nodes.dz(1) +1; % <-- viene decimale, corretto??
idL                     = sum(P.nodes.z < B.dL);

for i=1:P.kk
    tm                  = P.time(P.j) - B.tqstar(i);

%     CNH4_UR_CUM(i) = B.Cstar.UR(i)*(1-exp(-B.KhUR*tm));
%     CNH3_UR_CUM(i) = B.Cstar.UR(i)*(1-exp(-B.KvUR*tm));      

    
    % Calcolo del fattore di riduzione del coefficiente di nitrificazione
    % (vedi Gusman et al, 1999. Analytical Modeling of Nitrogen Dynamics in
    % soil and groundwater. J.of Irrigation and Drainage Engineering.
    for k=1:idL
        if P.teta(k)<=P.sh.tetafc(1)
            Fmw         = Fmw + P.teta(k)/P.sh.tetafc(1);
        else
            Fmw         = Fmw + P.sh.tetafc(1)/P.teta(k);
        end
    end
    % Il fattore di riduzione non viene calcolato nodo per nodo ma viene
    % mediato sull'intero strato di interramento.
    Fmw                 = Fmw / idL;

    % PERCHE' SCRIVO IN B?? dovrei usare una variabile temporanea!!
    B.KmORG_rp          = 5.6*10^12 * exp(-9800/(B.Tstar(P.kk)+273))*Fmw;
    B.KmORG_sw          = 4.0*10^9  * exp(-8400/(B.Tstar(P.kk)+273))*Fmw;

    % qui applica il decadimento (asintotico, quanto viene prodotto):
    CNH4_UR             = CNH4_UR       + B.Cstar.UR(i)*B.KhUR*exp(-B.KhUR*tm);
    % NH3 contribuisce per lo più nel momento della
    % somministrazione (per cui dopo non fai la sum):
    CNH3_UR             = CNH3_UR       + B.Cstar.UR(i)*B.KvUR*exp(-B.KvUR*tm);
    CNH4_ORG_rp         = CNH4_ORG_rp   + B.Cstar.ORG.rp(i)*B.KmORG_rp*exp(-B.KmORG_rp*tm);
    CNH4_ORG_sw         = CNH4_ORG_sw   + B.Cstar.ORG.sw(i)*B.KmORG_sw*exp(-B.KmORG_sw*tm);
end

% Non si moltiplica per W.dt perchè dalla deirvata calcolata al for
% precedente si ottiene il CNH4 prodotto per unità di tempo che è quello
% che deve entrare nella equazione ADE.

% L'apporto della forma solida NH4 ed NO3 dura per l'intero periodo
% P.kk:P.kk+1. essendo in g/cm2/d, l'input è già nella forma richiesta
% dall'ADE.
% Questo viene poi moltiplicato per W.dt nella stessa equazione.
CNH4_SD                 = B.Cstar.NH.SD(P.kk);
CNO_SD                  = B.Cstar.NO.SD(P.kk);
% Costruzione dei pool:
POOL_NH4_SD             = CNH4_UR + CNH4_ORG_rp + CNH4_ORG_sw + CNH4_SD;
POOL_NO_SD              = CNO_SD; % [g cm-2]
% Una volta costituito, il POOL si assume distribuito per l'intero spessore
% di suolo B.dL ==> lo si divide per B.dL e si ottiene una
% concnetrazione in g/cm3 di suolo.
for i=1:P.nz
    if i<=idL
        % se sono nello spessore di interramento:
        CNH4_pn(i)      = POOL_NH4_SD / B.dL;
        CNO_pn(i)       = POOL_NO_SD  / B.dL;
    else
        % se sono sotto lo spessore di interramento:
        CNH4_pn(i)      = 0;
        CNO_pn(i)       = 0;
    end
end
% -----------------------------------------------------------------------
% mineralizzazione sostanza organica :: END
% -----------------------------------------------------------------------

%% ADE
% -----------------------------------------------------------------------
% soluzione equazione convezione dispersione :: START
% -----------------------------------------------------------------------
for sl=1:2
    jj                  = 1;
    fluxin              = fluxsurf;
    fluxout             = P.flux(2);
    C1top               = P.Cinput(sl);
    C1bot               = C1(2,sl);
%     lambdatop           = S.lambda(1);
    lambdaup            = (S.lambda(1)+S.lambda(1))/2; % = S.lambda(1);
    lambdalow           = (S.lambda(1)+S.lambda(2))/2;
    
    while jj<=W.nlay
        for i=P.nodes.cumsum(jj)+1:P.nodes.cumsum(jj+1)
            % Antonio: we define all X in S.X in input file, why we
            % need to set a value here (only at first node of each layer
            % exept the first one).
            % ...there is something not good in these statements...
            if i == P.nodes.cumsum(jj+1)
%                 C1(i+1,sl)          = C1(i,sl); % non corretto farlo per ogni layer, ma solo con l'ultimo !!
                % CHECK if these 4 are valid and must be here!!
%                 S.lambda(i+1)   = S.lambda(i); % USED but I modified at the bottom of file. 
%                 S.Knitr(i+1)    = S.Knitr(i);  % UNUSED at i+1
%                 S.Kimmob(i+1)   = S.Kimmob(i); % UNUSED at all
%                 S.Kdenitr(i+1)  = S.Kdenitr(i);% UNUSED at all
            end
            
            % NEW VERSION
            fluxup              = ((P.flux(i)+fluxin )/2);
            fluxlow             = ((P.flux(i)+fluxout)/2);
            C1up                = (C1(i,sl)+C1top)/2;
            C1low               = (C1(i,sl)+C1bot)/2;
            % error: you are writing at P.nz+1 --> impossible!!
%             P.flux(P.nz+1)          = O.fluxbot(1,P.j,mm); % NOT GOOD!! --> DELETE!!
%             lambdaup                = (S.lambda(1)+lambdatop)/2;
%             lambdalow               = (S.lambda(1)+S.lambda(2))/2;
%             if i==P.nodes.cumsum(jj)+1      % FIRST NODE of LAYER jj
%                 fluxup              = ((P.flux(i)+fluxin)/2);
%                 fluxlow             = ((P.flux(i)+P.flux(i+1))/2);
%                 C1up                = (C1(i,sl)+C1top)/2;         
%                 C1low               = (C1(i,sl)+C1(i+1,sl))/2;
%                 if i>1                      % FIRST NODE of LAYERs > 1
%                     lambdaup        = (S.lambda(i)+S.lambda(i-1))/2;                
%                     lambdalow       = (S.lambda(i)+S.lambda(i+1))/2;
%                 end
%             else                            % INTERMEDIATE NODES of LAYER
%                 fluxup              = ((P.flux(i)+P.flux(i-1))/2);
%                 fluxlow             = ((P.flux(i)+P.flux(i+1))/2);
%                 C1up                = (C1(i,sl)+C1(i-1,sl))/2;                
%                 C1low               = (C1(i,sl)+C1(i+1,sl))/2;
%                 lambdaup            = (S.lambda(i) + S.lambda(i-1))/2;                
%                 lambdalow           = (S.lambda(i) + S.lambda(i+1))/2;
%             end

            % The following is the discretised Eq. for Cphys (tracer):
            % ...
            C2one                   = (fluxup*C1up - fluxlow*C1low)/P.nodes.dz(i);
            if i==P.nodes.cumsum(jj)+1
                C2two               = lambdaup*abs(fluxup)*(C1top-C1(i,sl)) / P.nodes.dz(1);
            else                
%                 C2two               = lambdaup*abs(fluxup)*(C1(i-1,sl)-C1(i,sl)) / P.nodes.dz(i);
                C2two               = lambdaup*abs(fluxup)*(C1top-C1(i,sl)) / P.nodes.dz(i);
            end
%             C2three                 = lambdalow*abs(fluxlow)*(C1(i,sl)-C1(i+1,sl))/P.nodes.dz(i);
            C2three                 = lambdalow*abs(fluxlow)*(C1(i,sl)-C1bot)/P.nodes.dz(i);
            C2_phys                 = W.dt*(-C2one + ((C2two-C2three)/P.nodes.dz(i)));
            % [g cm-3 H2O]
            C2phys                  = ( C2_phys + C1(i,sl)*P.teta(i) ) / P.teta(i);
            
            if P.teta(i)<P.sh.tetafc(i)
                teta_ratio          = P.teta(i)/P.sh.tetafc(i);
            else
                teta_ratio          = P.sh.tetafc(i)/P.teta(i);
            end
            C2_ntf_lq               = S.Knitr(i)*1.07^ ...
              (B.Tstar(P.kk)-S.Topt)*teta_ratio*C1(i,1)*P.teta(i);
            C2_ntf_sd               = S.Knitr(i)*1.07^ ...
              (B.Tstar(P.kk)-S.Topt)*teta_ratio*P.sh.dap(i)*S1(i,1);
            % attingimento selettivo:
            C2_sink                 = S.NX.Kr(sl)*P.sink(i)*C1(i,sl);

            if sl==1
                % SsSk: Source-Sink for solutes [g cm-3 suolo day-1]
                SsSk                = -C2_ntf_lq-C2_ntf_sd-C2_sink+CNH4_pn(i);
            elseif sl==2
                % Sm(z,t), Eq. 16
                C2_immb             = S.Kimmob(i)*1.05^ ...
                  (B.Tstar(P.kk)-S.Topt)*teta_ratio*C1(i,2)*P.teta(i);
                
                % teta_tsh: the threshold water content for de-nitrification:
                teta_tsh            = 0.627*P.sh.tetafc(i)-0.0267*(P.sh.tetas(i) ...
                                     -P.teta(i))/P.sh.tetas(i)*P.sh.tetafc(i);
                if P.teta(i)>teta_tsh
                    % Sd(z,t), Eq. 17
                    C2_dntf         = S.Kdenitr(i)*1.07^(B.Tstar(P.kk)- ...
                                      S.Topt)*(P.teta(i)-teta_tsh) / ...
                                      (P.sh.tetafc(i)-teta_tsh)*C1(i,2)*P.teta(i);
                else
                    C2_dntf         = 0;
                end
                % SsSk: Source-Sink for solutes [g cm-3 suolo day-1]
                SsSk                = +C2_ntf_lq +C2_ntf_sd -C2_immb ...
                                      -C2_dntf -C2_sink +CNO_pn(i);
            end
            
            % è giusto aggiungere C1?? Risponde Antonio...
            C2_chem                 = ( SsSk*W.dt + C1(i,sl)*P.teta(i) ) ...
                                      / P.teta(i);
            C2_tot                  = C2phys + C2_chem;

            if S.NX.Kf2(sl)==1
                C2out(i,sl)   = ( C2phys*P.teta(i) + P.sh.dap(i)* ...
                                      S.NX.Kf1(sl)*C1(i,sl) + SsSk*W.dt ) ...
                                      / ( fnteta( W.hfc, P.sh, i ) + ...
                                      P.sh.dap(i)*S.NX.Kf1(sl) );
            elseif C2_tot < 10^-9
                C2phys              = 0;
                C2out(i,sl)   = C2phys;    
            else
                % **Bisection method**
                % Bisection method (vedi libro con esercizi numerici in
                % Matlab) for adsorbing solutes with Freundlich f_bis_sol è
                % la funzione di cui si vuole trovare lo zero.
                % E' bene fissare un limite inferiore per
                % C2out(i,sl)=10^k per eliminare eventuali valori
                % negativi ed evitare problemi di convergenza del metodo.
                f_bis_sol           = @(yps) ...
                        P.teta(i)*yps        - ...
                        P.teta(i)*C2phys     + ...
                        P.sh.dap(i)*S.NX.Kf1(sl)*yps^S.NX.Kf2(sl) - ...
                        P.sh.dap(i)*S1(i,sl)          - ...
                        SsSk*W.dt;
                % bisection calc:
                C2out(i,sl)   = mln_bisection( f_bis_sol, -100, 100, -30 );
            end
            %--------------------------------------------------------------
            % ***NEW VERSION***
            % operazioni prima di passare al prossimo strato:
            % (remember that this is the i node and I'll manage the i+1
            %  node)
            fluxin                  = P.flux(i+0);          % i-1
            switch i
                case P.nz
                    % do nothing! I'm exiting the loop
                case P.nz-1
                    fluxout         = fluxbot;
                    C1top           = C2out(i+0,sl);        % i-1
                    C1bot           = C1(i+1,sl);           % i+0        **particular case!
                    lambdaup        = (S.lambda(i+1)+S.lambda(i+0))/2;
                    lambdalow       = (S.lambda(i+1)+S.lambda(i+1))/2; % **particular case!
                case P.nz-2
                    fluxout         = fluxbot;
                    C1top           = C2out(i+0,sl);        % i-1
                    C1bot           = C1(i+1,sl);           % i+0        **particular case!
                    lambdaup        = (S.lambda(i+1)+S.lambda(i+0))/2;
                    lambdalow       = (S.lambda(i+1)+S.lambda(i+1))/2; % **particular case!                    
                otherwise
                    fluxout         = P.flux(i+2);          % i+1
                    C1top           = C1(i+0,sl);           % i-1
                    C1bot           = C1(i+2,sl);           % i+1
                    lambdaup        = (S.lambda(i+1)+S.lambda(i+0))/2;
                    lambdalow       = (S.lambda(i+1)+S.lambda(i+2))/2;
            end
            %--------------------------------------------------------------
        end
        
        %------------------------------------------------------------------
        % ***OLD VERSION***
        % operazioni prima di passare al prossimo strato:
        % CORRECT IT !!
        % error: you are reading at (jj+1)+1 --> it should be (jj+1),
        % infact you use P.flux(i-1) as fluxin for current node i above (at
        % the beginning of for i=...)
        %   fluxin        = P.flux(P.nodes.cumsum(jj+1)+1);
        %   C1top         = O.C2(  P.nodes.cumsum(jj+1)+1,P.j,mm,sl);
        % CORRECT IT !!
        
        jj                                  = jj+1;
        %-----------------------------------------------------------------
    end
end
% -----------------------------------------------------------------------
% soluzione equazione convezione dispersione :: END
% -----------------------------------------------------------------------
return