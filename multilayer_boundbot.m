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
switch W.botbc
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
%         error('Bottom boundary condition (W.botbc=%d) not yet implemented',W.botbc)
        dFdh_botB = 0;
    
    case 8% 8  Free outflow at soil-air interface
        error('Bottom boundary condition (W.botbc=%d) not yet implemented',W.botbc)
    
    otherwise
        error('Bottom boundary condition (W.botbc=%d) not properly set',W.botbc)
end