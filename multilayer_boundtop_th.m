%% check fluxes
%             P.tq                = B.top.thqstar(P.kk+1);
%% ---top       boundary
% if W.itbc==0                    % flux
%     W.qsurf         = B.top.hqstar(P.tidx);
% elseif and(W.itbc==1,P.rnf==0)  % potential
%     W.hsurf         = B.top.hqstar(P.tidx); 
% end

% restore qsurf & hsurf -- check with Antonio!! --> put in boundtop_th.m!!
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
%   OR
% to change top boundary condition on a new dt (when no convergence was
% reached)
if W.itbc==0
    W.qsurf     = B.top.hqstar(P.tidx);
elseif and(W.itbc==1,P.rnf==0)
    W.hsurf     = B.top.hqstar(P.tidx);
elseif W.itbc==1 && P.rnf==1 && P.L==1
    W.qsurf     = B.top.hqstar(P.tidx);
end
%% ---bottom    boundary
if W.ibbc==0                % flux
    W.qbot          = P.bothq(P.tidx);
elseif W.ibbc==1            % potentials
    W.hbot          = P.bothq(P.tidx);
end
%% solutes
if and(W.isol==2,W.iCtopvar==1)% check with Antonio!!
    % fertigation(NH4):
    P.Cinput(1)     = P.compounds(6,P.tidx);
    % fertigation(NO3):
    P.Cinput(2)     = P.compounds(7,P.tidx);
end
%% crop
if W.iveg==1
%     P.ETp0          = P.Kc(P.tidx)*B.top.eto(P.tidx);
    P.ETp0          = P.Kc(P.tidx)*V.ETr(P.tidx);
    P.Ep            = P.ETp0*exp(-V.extf*V.LAI(P.tidx));
    P.Tp            = P.ETp0-P.Ep;
    P.Droot         = V.Droot(P.tidx);
end