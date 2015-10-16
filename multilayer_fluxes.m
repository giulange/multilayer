%% ?? qimmob ??
% qimmob:       Soil water flux between mobile and immobile fraction in
%               case of fingered flow [cm d-1].
qimmob      = zeros(P.nz,1);% unable to understand where it's defined in SWAP-32!

% QExcMpMtx:    Water exchange flux between matrix and macropores per
%               compartment [cm d-1].
QExcMpMtx   = zeros(P.nz,1);% defined in MacroRate.for:71,294
%% Flag indicating precribed groundwater level below bottom soil column:
%   -not implemented yet;
fllowgwl    = 0;
%% TO BE ACTIVATED:
% --- determine qbot if not specified
% if W.swbotb==5 || W.swbotb==7 || W.swbotb==8 || W.swbotb==-2 || ...
%         (W.swbotb==1 && fllowgwl)
%     qbot = qtop + qrosum + qdrtot - QMaPo + (volact-volm1)/P.dt;
% end
%% ADAPTED from SWAP-32

% --- calculate fluxes [cm d-1] from changes in volume per compartment
P.q(P.nz+1) = P.qbot;
% inq(i) = inq(i) + q(i)*P.dt;
for ii = P.nz:-1:1
    P.q(ii) = -(P.teta(ii)-P.teta_jm1(ii) + qimmob(ii)) ...
              * P.FrArMtrx(ii)*P.nodes.dz(ii)/P.dt + ...
              + P.q(ii+1)-P.sink(ii)+QExcMpMtx(ii);
% I used qrot = P.sink above!

%     for level=1:nrlevs
%        q(ii) = q(ii) - qdra(level,ii);
%     end
%     inq(i) = inq(i) + q(i)*P.dt;
end
%% DARCY ~ A.Coppola
% P.q(P.nz)                = -P.K(P.nz) * ((P.h(P.nz-1)-W.hbot) ...
%                               /(0.5*P.nodes.dz(P.nz+1)+P.nodes.dz(P.nz))+1);
% 
% P.q(2:P.nz-1)            = -P.K(2:P.nz-1).*( ...
%                               ( P.h((1:P.nz-2))-P.h((3:P.nz)) ) ...
%                               ./ (P.nodes.dz(2:P.nz-1)+0.5*P.nodes.dz(1:P.nz-2)+0.5*P.nodes.dz(3:P.nz))+1 );
% if W.itbc==1 || W.qsurf==P.Emax
%     P.q(1)               = -P.K(1) * ((W.hsurf-P.h(2)) /(1.5*P.nodes.dz(1))+1);
% else
%     % check this code!!
%     P.q(1)               = (O.fluxsurf(1,P.j,mm)*(1/(0.5*P.nodes.dz(1))) + ...
%                               P.q(2)*(1/(0.5*(P.nodes.dz(1)+P.nodes.dz(2)))) ...
%                               ) / ( (1/(0.5*P.nodes.dz(1)))+ 1/(0.5*(P.nodes.dz(1)+P.nodes.dz(2))) );
% end