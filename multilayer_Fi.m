%% NOTES
% This script calculates the function F used to solve the Newton-Raphson
% p-iteration.
% **hsurf/qtop: develop a multilayer_boundary_top.m to manage top boundary:
% *** I need to define Fi at top and bottom nodes !!!!
% see headcalc.for:
%   -F(1) at line 442:
%       F(1) = (theta(1)-thetm1(1))*FrArMtrx(1)*dz(1)/dt + sink(1) + qrot(1) + kmean(2) * hgrad(2)
%       if(ftoph)then
%       	hgrad(1) = (hsurf-h(1))/disnod(1) + 1.d0
%           F(1) = F(1) - kmean(1) * hgrad(1)
%       else
%           F(1) = F(1) + qtop
%       end if
%    where:
%       FrArMtrx    --> Fraction of horizontal area of soil matrix per
%                       compartment [-]
%       ftoph       --> Flag indicating that the pressure head is
%                       prescribed at the soil surface
%       kmean(2)    --> Should be Kip2(2)
%       hgrad(2)    --> Should be (h_p(1)-h_p(2))./(dz(1)+dz(2))
%% init
Fi          = zeros(P.nz,1);
%% bottom bounday
% **botB: develop a multilayer_boundary_bottom.m to manage bottom boundary:
%       
%	-F(NN) at line 493
%  
botB = 0;
%% ** Fi **
% Function Fi (SWAP 32 manual, Eq. 2.28, page 36)
%   -It computes Fi before entering Newton-Raphson iteration;
%   -check if its computation is different from that within p-iteration;
%% *top          node
% [headcalc.for, lines 175-180]
% Determines whether the top boundary is flux or head controlled:
%   -OLD :: W.itbc (1-->(~isflux), 0-->isflux)
if ~isflux%     HEAD controlled top boundary
    % Fi-function:
    Fi_topB     = - 2.*Kim2(1).*(hsurf-h_pm1(1))./(dz(1)) - Kim2(1);
elseif isflux%  FLUX controlled top boundary
    % Fi-function:
    Fi_topB     = qtop;
end
% [headcalc.for, lines 161-162]
Fi(1)   = dz(1)./P.dt.*(teta_p(1) - teta_jm1(1)) + sink(1) + macr(1) + ...
        + Fi_topB + ...
        + 2.*Kip2(1).*(h_p(1)-h_p(2))./(dz(1)+dz(2)) + Kip2(1) ;
%% *bottom       node
% *TO BE IMPLEMENTED:
Fi_botB = 0;
warning('develop the value of Fi_botB!!')
Fi(P.nz)= dz(1)./P.dt.*(teta_p(P.nz) - teta_jm1(P.nz)) + sink(P.nz) + macr(P.nz) + ...
        - 2.*Kim2(P.nz).*(h_p(P.nz-1)-h_p(P.nz))./(dz(P.nz-1)+dz(P.nz)) - Kim2(P.nz) + ...
        + Fi_botB;
warning('Fi_botB not implemented yet!!')
%% *intermediate nodes
Fi(i_i) = dz_i./P.dt.*(teta_p(i_i) - teta_jm1(i_i)) + sink(i_i) + macr(i_i) + ...
        - 2.*Kim2(i_i).*(h_p(i_im1)-h_pm1(i_i))./(dz_im1+dz_i) - Kim2(i_i) + ...
        + 2.*Kip2(i_i).*(h_p(i_i)-h_p(i_ip1))./(dz_i+dz_ip1) + Kip2(i_i) ;
%