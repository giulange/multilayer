% --- determine qbot if not specified
%% Flag indicating precribed groundwater level below bottom soil column:
%   -not implemented yet;
fllowgwl = 0;

%% TO BE ACTIVATED:
% if W.swbotb==5 || W.swbotb==7 || W.swbotb==8 || W.swbotb==-2 || ...
%         (W.swbotb==1 && fllowgwl)
%     qbot = qtop + qrosum + qdrtot - QMaPo + (volact-volm1)/P.dt;
% end

%% ORIGINAL SWAP-32
% % --- calculate fluxes [cm d-1] from changes in volume per compartment
% P.q(P.nz+1) = P.qbot;
% % inq(i) = inq(i) + q(i)*P.dt;
% for ii = P.nz:-1:1
%     P.q(ii) = -(P.teta(ii)-P.teta_jm1(ii) + ...
%               + qimmob(ii))*P.FrArMtrx(ii)*P.nodes.dz(ii)/P.dt + ...
%               + P.q(ii+1)-qrot(ii)+QExcMpMtx(ii);
% 
%        
% %     for level=1:nrlevs
% %        q(ii) = q(ii) - qdra(level,ii);
% %     end
% %     inq(i) = inq(i) + q(i)*P.dt;
% 
% end
%% DARCY
P.flux(P.nz)                = -P.K(P.nz) * ((P.h(P.nz-1)-W.hbot) ...
                              /(0.5*P.nodes.dz(P.nz+1)+P.nodes.dz(P.nz))+1);

P.flux(2:P.nz-1)            = -P.K(2:P.nz-1).*( ...
                              ( P.h((1:P.nz-2))-P.h((3:P.nz)) ) ...
                              ./ (P.nodes.dz(2:P.nz-1)+0.5*P.nodes.dz(1:P.nz-2)+0.5*P.nodes.dz(3:P.nz))+1 );
if W.itbc==1 || W.qsurf==P.Emax
    P.flux(1)               = -P.kond(1) * ((W.hsurf-P.h1(2)) /(1.5*P.nodes.dz(1))+1);
else
    % check this code!!
    P.flux(1)               = (O.fluxsurf(1,P.j,mm)*(1/(0.5*P.nodes.dz(1))) + ...
                              P.flux(2)*(1/(0.5*(P.nodes.dz(1)+P.nodes.dz(2)))) ...
                              ) / ( (1/(0.5*P.nodes.dz(1)))+ 1/(0.5*(P.nodes.dz(1)+P.nodes.dz(2))) );
end
