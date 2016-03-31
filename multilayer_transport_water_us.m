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
%% init
% Newton-Raphson/Thomas was broken by an exit condition:
nr_breaked  = false;
% Newton-Raphson/Thomas convergence not possible at this timestep:
fl_noconv   = false;
% Newton-Raphson/Thomas number of non-convergences producing dt reduction:
n_noconv    = 0;
%% Module for water transport
switch W.wt_mod
    case 0
%% *THOMAS
        % run only if P.timestep raise of 1 integer step!
%         if P.L==1
%             multilayer_boundtop_th
%         end
        while ~nr_breaked && ~fl_noconv
            % *TOP
            multilayer_boundtop_th
            
            % *THOMAS
%             run multilayer_wtus_thomas_OLD_fn.m
            multilayer_wtus_thomas_OLD_fn
            % [P.iter(:,P.j-5:P.j);P.time(P.j-5:P.j)-P.time(P.j-6:P.j-1)]
        end
        P.iter(1:4,P.j)     = [p;nr_breaked;fl_noconv;n_noconv];
        % *RUNOFF [& runon?]
        O.runoff(1,P.j,mm)  = W.qsurf-O.fluxsurf(1,P.j,mm);
        P.h = O.h22(:,P.j,mm);
        for i=1:P.nz
            P.teta(i) = fnteta_OLD(P.h(i), P.sh.tetas(i),P.sh.tetar(i),P.sh.alfrs(i),P.sh.fi(i),P.sh.alfvg(i),P.sh.en(i),P.sh.alfvg2(i),P.sh.en2(i),P.sh.ifr(i));
        end
        % *W.itbc
        %   Questo controllo va fatto solo se W.itopvar=1.
        %   Se W.itopvar=0, una volta avvenuto il cambio da W.itbc=0 a
        %   W.itbc=1, non e' piu' necessario tornare ad W.itbc=0
        %   perche' il flusso non cambia fino a fine simulazione.
        %   Se il nuovo B.top.hqstar(P.kk) e' maggiore del vecchio
        %   (B.top.hqstar(P.kk-1) non e' necessario cambiare W.itbc da
        %   1 a 0.
        %   Se invece e' minore allora bisogna di nuovo verificare se
        %   il nuovo B.top.hqstar sia maggiore di P.fluxsurf_max, ed
        %   allora occorre cambiare W.itbc da 1 a 0.
        %   Ovviamente questa verifica va fatta solo se si è entrati in
        %   W.itbc=1 partendo da W.itbc=0.
        if P.L==1               % flag counter top-bound
            % se il flusso al nuovo P.tidx > P.tidx-1 ho ancora runoff
            if and(P.rnf==1,abs(B.top.hqstar(P.tidx))<abs(B.top.hqstar(P.tidx_jm1)))
                P.rnf=0;        % potrebbe non avere piu' senso 
                W.itbc=0;
            end
        end
    case 1
%% *NEWTON-RAPHSON
        % bottom boundary condition:
        multilayer_boundbot

        % SAVE:
        %   - This SAVE must be put inside or outside the while? --> SWAP
        %     put it outside the loop!
        %       (1) If it is inside I hold the h from previous
        %           iteration step to start a new iteration based on a
        %           smaller dt (when no convergence was reached).
        %       (2) If it must be put outside, then I have to update
        %           h, teta and pond at both j-1 and j to the values
        %           stored at previous j (i.e. h=O.h22(:,P.j-1,mm) and
        %           so on).
        % Pressure head at former time level [cm]
        P.h_jm1 = P.h;
        % Volumic soil water content at former time level [-]
        P.teta_jm1 = P.teta;
        % GroundWater Level:
        %gwl_jm1 = gwl;
        % Ponding:
        pond_jm1 = pond;

        % Richard's equation solution by means of Newton-Raphson:
        while ~nr_breaked && ~fl_noconv            
            multilayer_wtus_newtonraphson
        end

        % STORE:
        P.iter(:,P.j)           = [p;nr_breaked;fl_noconv;n_noconv;iL;bt_breaked];
        O.h22(:,P.j,mm)         = P.h;%     *store current PRESSURE HEAD [cm]
        O.runoff(1,P.j,mm)      = runoff;%  *store current RUNOFF   [cm]
        O.fluxsurf(1,P.j,mm)    = P.qtop;%  *store current QTOP     [cm d-1]
        O.pond(1,P.j,mm)        = pond;%    *store current PONDING  [cm]
        O.fluxbot(1,P.j,mm)     = P.qbot;%	*store current QBOT     [cm d-1]
        % REGISTER: {K,K1s2,flux,macropore}
        %   -same as SWAP-32 call to SoilWater(3)
        multilayer_soilwater_register
end
%% clean
clear p nr_breaked fl_noconv n_noconv iL bt_breaked