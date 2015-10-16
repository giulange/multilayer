%% ** Fi **
% This script calculates the function F used to solve the Newton-Raphson
% p-iteration.
% 	-SWAP 32 manual, Eq. 2.28, page 36
%   -include k for implicit/explicit
%% hgrad
%   -in SWAP-32 ¶disnod(i) = 0.5*(dz(i-1)+dz(i))¶
hgrad           = NaN(P.nz,1);
hgrad(2:P.nz)   = (P.h(1:P.nz-1) - P.h(2:P.nz)) ./ P.nodes.disnod(2:P.nz) +1;
%% *top          node
% [headcalc.for, lines 175-180]
if ~isflux%     HEAD controlled top boundary
    % Fi-function:
    Fi_topB     = - P.Kim2(1).*((hsurf-P.h(1))./P.nodes.disnod(1) +1 );
elseif isflux%  FLUX controlled top boundary
    % Fi-function:
    Fi_topB     = P.qtop;
end
% [headcalc.for, lines 161-162]
Fi(1)   = P.nodes.dz(1)./P.dt.*(P.teta(1) - P.teta_jm1(1))*P.FrArMtrx(1) + ...
        + P.sink(1) + macr(1) + ...
        + Fi_topB + ...
        + P.Kip2(1) .* hgrad(2);
%% *intermediate nodes
Fi(2:P.nz-1) = P.nodes.dz(2:P.nz-1)./P.dt.*(P.teta(2:P.nz-1) + ...
             - P.teta_jm1(2:P.nz-1)).*P.FrArMtrx(2:P.nz-1) + ...
             + P.sink(2:P.nz-1) + macr(2:P.nz-1) + ...
             - P.Kim2(2:P.nz-1).*hgrad(2:P.nz-1) + ...
             + P.Kip2(2:P.nz-1).*hgrad(3:P.nz);
%% *bottom       node
if P.j==1
%     warning('%s (%s)!\n  %s\n\n','Delete the following switch (W.SwBotB) from within THIS script', ...
%             mfilename, 'It is already computed in multilayer_boundtop.m !!');
end
switch W.SwBotB
    case 1% 1  Prescribe groundwater level
    case 2% 2  Prescribe bottom flux
    case 3% 3  Calculate bottom flux from hydraulic head of deep aquifer
    case 4% 4  Calculate bottom flux as function of groundwater level
    case 5% 5  Prescribe soil water pressure head of bottom compartment
    case 6% 6  Bottom flux equals zero
    case 7% 7  Free drainage of soil profile
        if W.SwMacro
            P.Kip2(P.nz) = P.Kip2(P.nz)*P.FrArMtrx(P.nz);
        end
        P.qbot      = -1 * P.Kip2(P.nz);
        Fi_botB     = -P.qbot;

    case 8% 8  Free outflow at soil-air interface
    otherwise
        error('Bottom boundary condition (W.botbc=%d) not properly set',W.botbc)
end

Fi(P.nz)= P.nodes.dz(P.nz)./P.dt.*(P.teta(P.nz) - P.teta_jm1(P.nz))*P.FrArMtrx(P.nz) + ...
        + P.nodes.dz(P.nz) * (P.sink(P.nz) + macr(P.nz)) + ...
        - P.Kim2(P.nz).*hgrad(P.nz) + ...
        + Fi_botB;
%
%% *macropore
if W.SwMacro
    Fi = Fi - QExcMpMtx;
end