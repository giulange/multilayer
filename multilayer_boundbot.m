%% OPTIONS that can be implemented
% 1  Prescribe groundwater level
% 2  Prescribe bottom flux
% 3  Calculate bottom flux from hydraulic head of deep aquifer
% 4  Calculate bottom flux as function of groundwater level
% 5  Prescribe soil water pressure head of bottom compartment
% 6  Bottom flux equals zero
% 7  Free drainage of soil profile
% 8  Free outflow at soil-air interface
%%
switch W.SwBotB
    case 1% 1  Prescribe groundwater level
         error('Bottom boundary condition (W.botbc=%d) not yet implemented',W.botbc)
    case 2% 2  Prescribe bottom flux
        error('Bottom boundary condition (W.botbc=%d) not yet implemented',W.botbc)
    case 3% 3  Calculate bottom flux from hydraulic head of deep aquifer
        error('Bottom boundary condition (W.botbc=%d) not yet implemented',W.botbc)
    case 4% 4  Calculate bottom flux as function of groundwater level
        error('Bottom boundary condition (W.botbc=%d) not yet implemented',W.botbc)
    case 5% 5  Prescribe soil water pressure head of bottom compartment
        error('Bottom boundary condition (W.botbc=%d) not yet implemented',W.botbc)
    case 6% 6  Bottom flux equals zero
        error('Bottom boundary condition (W.botbc=%d) not yet implemented',W.botbc)
    case 7% 7  Free drainage of soil profile
        % This case is called as unitary-gradient by A.Coppola, who applies
        %   qbot = -K * ( (h(nz)-hbot)/dz(nz+1) +1 )
        % where
        %   hbot = h(nz) - (W.grad-1)*dz(nz+1)
        % with W.grad=1 is:
        W.hbot      = P.h(P.nz);
        % which, for W.grad=1, reduces to the qbot computed below:
        % *Fi
        P.qbot      = -1 * P.Kip2(P.nz);
%         figure(1),pause(0.1),plot( P.q(101-5:101),1:6 )
    case 8% 8  Free outflow at soil-air interface
        error('Bottom boundary condition (W.botbc=%d) not yet implemented',W.botbc)
    otherwise
        error('Bottom boundary condition (W.botbc=%d) not properly set',W.botbc)
end