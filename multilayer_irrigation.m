%% NOTES
% 1.\   When computing irrigation volume we should also account for the
%       irrigation type/system
% 2.\   How to include salinity/concentration (cirr) of irrigation water?
% 3.\   How to include depth (rDEPTH) at which irrigation water is given?
%% A.Coppola
%%controllo irrigazione a domanda (irr_dem: 0=no irrigation on demand;%%1=irrigation on demand)
%%controlla che il potenziale medio nella root zone sia inferiore al potenziale critico per l'irrigazione (hcrit).
%%una volta raggiunto l'hcrit, calcola la differenza di storage (delta_stor) fra la capacità di campo ed il contenuto 
%%d'acqua attuale.Quindi calcola il valore del potenziale flusso irriguo dividendo delta_stor per il deltat fra il tempo 
%%attuale (tqstra(kk)e quello del nuovo input (tq). Se qirr(kk) è maggiore del qtstar(kk), la differenza (delta_qirr(kk)) 
%%viene aggiunta al qtstar(kk) aggiornandolo anche per i tempi successivi fino al prossimo kk. Il calcolo di deltaqirr 
%%sembra superfluo ma serve a conservare memoria del volume irriguo somministrato
%%il controllo va fatto solo nei periodi in cui il suolo è effettivamente 
%%coltivato. quando il suolo è nudo il controllo va escluso

% if irr_dem==1
%%yc=5 years block counter (5 years=1825 days. In pratica i periodi nei
%%quali il suolo è nudo si ripetono negli stessi intervalli ogni 5 anni (nel caso studio considerato)
%%per esempio, nel primo blocco di 5 anni gli intervalli sono 59-89,244-271 e 456-819. nel secondo blocco di anni si ripetono 
%%negli stessi intervalli ma traslati nel tempo di 5 anni (1825 days) 
%     if time(j)>1825*yc; 
%         yc=yc+1;
%     end
%     %%intervalli temporali suolo nudo
%     itsn = [59:89,244:271,456:819];
%     
%     if isempty (intersect( floor(mod( time(j), 1825*(yc-1) )) , itsn ))
%         if and(itbc==0,T==1)
%             hrz=mean(h1(7:22,j));
%             if hrz<hcrit
%                 for i=1:22
%                     dstor(i)=(tetafc(i)-teta(i,j))*dz(i);
%                 end
%                 delta_stor=sum(dstor);
%                 qirr(kk)=-delta_stor/(tq-tqstar(kk)); 
%                 if qirr(kk)<=qtstar(kk);
%                     %                       hrz, qtstar(kk)
%                     delta_qirr(kk)=qirr(kk)-qtstar(kk);    
%                     qtstar(kk)=qtstar(kk)+delta_qirr(kk);
%                     qtstar(kk)
%                 end
%             end
%         end
%     end
% end
%% VOLUME
switch V.flirri
    case 0
        error('You should not be here with V.flirri=0')
        
    case 1% schedule:known,   volume:known
    % **THOMAS ALGORITHM
        % B.top.hqstar is updated in multilayer_boundtop_th
    % **NEWTON-RAPHSON ALGORITHM
        % What to do and where?
        error('You should not be here with V.flirri=1')
        
    case 2% schedule:known,   volume:unknown
        % insert block by A.Coppola, parameterised on h_to
        error('Not implemented yet!')
        
    case 3% schedule:unknown, volume:unknown
        % insert block by A.Coppola, parameterised on h_from & h_to
        if isnan( P.h_from(P.tidx) ), break, end
        
        % root zone usefull for irrigation:
        rz_irri             = P.Droot(P.tidx) * V.DrootE;
        % *TOP node
        Fnod = find(P.nodes.z <= rz_irri(1),1,'last');
        if isempty(Fnod)
            rz_node(1)      =       +1;
        else
            rz_node(1)      = Fnod  +1;
        end
        % *BOTTOM node
        Fnod = find(P.nodes.z <= rz_irri(2),1,'last');
        if isempty(Fnod)
            rz_node(2)      =       +1;
        else
            rz_node(2)      = Fnod  +1;
        end
        i                   = rz_node(1):rz_node(2);
        % **WATER STRESS:
        % Critical root zone pressure head, useful for irrigation requirement:
        P.hrz_cm(P.tidx)    = sum(O.h22(i,P.j-1) .* P.nodes.dz(i)) / sum( P.nodes.dz(i) );
        % if the pressure head threshold is overcame:
        if P.hrz_cm(P.tidx) <= P.h_from(P.tidx)
            i               = 1:rz_node(2);
            % The volumetric soil water content at prescribed h_to
            % threshold, and at nodes of soil grid interested by root zone
            % used to calculate crop stress:
            %teta_to         = fnteta( repmat(P.h_to(P.tidx),length(i),1), P.sh, i );
            teta_to         = multilayer_fnteta_vgm( repmat(P.h_to(P.tidx),length(i),1), P.sh, i );
            % **WATER REQUIREMENT:
            % amount of irrigation water needed to fill the soil till the
            % h_to threshold, expressed in [cm]:
            P.dstor(P.tidx) = (teta_to - P.teta(i))' * P.nodes.dz(i);
            % the net flux of irrigation:
            netqirri        = +P.dstor(P.tidx) / W.timestep -B.top.rain(P.tidx);
            % [cm]
            if netqirri > 0.5%cm
                % **WATER GIVEN:
                P.irri(P.tidx) = netqirri;
            end
        end
        
    otherwise
        error('Bad V.flirri(=%d) configuration',V.flirri)
end