%% NOTES
% This script computes the partial derivatives of Fi to the pressure head:
%   dF/dh = ...
% The following Jacobian terms are computed with reference to the SWAP-32
% implementation:
%  *MAIN diagonal (SWAP-32 has dFdhM)
%    top            dF_dh(1,1)
%    intermediate   dF_dh(2:nz-1,2:nz-1)
%    bottom         dF_dh(nz,nz)
%  *ABOVE diagonal (SWAP-32 has dFdhL)
%  *BELOW diagonal (SWAP-32 has dFdhU)
%% PRE
dz              = P.nodes.dz;
disnod(1)       = 0.5*dz(1);
disnod(2:P.nz,1)= 0.5*(dz(1:P.nz-1) + dz(2:P.nz));
hgrad(1)        = (hsurf-P.h(1))/disnod(1) +1;
hgrad(2:P.nz)   = (P.h(1:P.nz-1) - P.h(2:P.nz)) ./ disnod(2:P.nz) +1;
%% *Derivative of the internodal hydraulic conductivity to the pressure
%   -this are computed in SWAP-32 as the product dkdh*dkmean
% head (SWAP-32 manual, Appendix 5, page 239):
dKim2_dhi(2:P.nz-1)   = multilayer_dK_is2_dh( dKi_dhi(2:P.nz-1),   W.Kmeth, P.K(1:P.nz-2), P.K(2:P.nz-1), dz(1:P.nz-2), dz(2:P.nz-1), false );
dKip2_dhi(2:P.nz-1)   = multilayer_dK_is2_dh( dKi_dhi(2:P.nz-1),   W.Kmeth, P.K(3:P.nz-0), P.K(2:P.nz-1), dz(3:P.nz-0), dz(2:P.nz-1), false );
dKim2_dhim1(2:P.nz-1) = multilayer_dK_is2_dh( dKi_dhi(1:P.nz-2),   W.Kmeth, P.K(1:P.nz-2), P.K(2:P.nz-1), dz(1:P.nz-2), dz(2:P.nz-1), true  );
dKip2_dhip1(2:P.nz-1) = multilayer_dK_is2_dh( dKi_dhi(3:P.nz-0),   W.Kmeth, P.K(3:P.nz-0), P.K(2:P.nz-1), dz(3:P.nz-0), dz(2:P.nz-1), true  );
% top    node particular cases:     check with Antonio
dKip2_dhi(1)     = multilayer_dK_is2_dh( dKi_dhi(1), W.Kmeth, P.K(2), P.K(1), dz(2), dz(1), false );
dKip2_dhip1(1)   = multilayer_dK_is2_dh( dKi_dhi(2), W.Kmeth, P.K(2), P.K(1), dz(2), dz(1), true  );
% bottom node particular cases:     check with Antonio
dKim2_dhi(P.nz)  = multilayer_dK_is2_dh( dKi_dhi(P.nz),   W.Kmeth, P.K(P.nz-1), P.K(P.nz), dz(P.nz-1), dz(P.nz), false );
dKim2_dhim1(P.nz)= multilayer_dK_is2_dh( dKi_dhi(P.nz-1), W.Kmeth, P.K(P.nz-1), P.K(P.nz), dz(P.nz-1), dz(P.nz), true  );
dKip2_dhi(P.nz) = dKi_dhi(P.nz) /2;% [headcalc.for:326, I'm not sure it's my case!!] -- check with Antonio
% partial derivative of the macro-pore exchange to the pressure head:
if W.SwMacro % && ~fl_UnsatOk(3)
    multilayer_macropore% --> outputs{ dSINKmacr_dh, ... }
end
%% Jacobian coefficients of the tri-diagonal system of equations
% (SWAP-32 manual, Appendix 4):
%   -Partial derivatives of Fi to pressure heads parameterised upon k
%   to allow for treating the hydraulic conductivities implicitly (when
%   k=1)
%% *MAIN diagonal
%%   *top     boundary :: dF_dh(1,1)
%       -OLD :: W.itbc (1-->(~isflux), 0-->isflux)
dF_dh(1,1)  = dz(1)/P.dt*Ci(1)*P.FrArMtrx(1) + dz(1)*dSINKmacr_dh(1) + ...
            + P.Kip2(1)/disnod(2) + ...
            + P.Kim2(1)/disnod(1) * (~isflux) + ...% HEAD controlled top boundary
        k*( - dKi_dhi(1)*hgrad(1) * (~isflux) + ...% HEAD controlled top boundary;
            + dKip2_dhi(1)*hgrad(2)...
          );%                                        k implicit/explicit
%%   *intermediate nodes
%     [dFi_dhi(2:P.nz-1,2:P.nz-1)]
IDX         = sub2ind([P.nz,P.nz],2:P.nz-1,2:P.nz-1);
dF_dh(IDX)  = dz(2:P.nz-1)/P.dt.*Ci(2:P.nz-1).*P.FrArMtrx(2:P.nz-1) + ...
            + dz(2:P.nz-1).*dSINKmacr_dh(2:P.nz-1)  + ...
            + P.Kim2(2:P.nz-1)./disnod(2:P.nz-1)    + ...
            + P.Kip2(2:P.nz-1)./disnod(3:P.nz-0)    + ...
        k*( - dKim2_dhi(2:P.nz-1).*hgrad(2:P.nz-1)  + ...
            + dKip2_dhi(2:P.nz-1).*hgrad(3:P.nz-0)    ...
          );
%%   *bottom  boundary:
switch W.SwBotB
    case {1,2,3,4,5,6}
        error('not implemented yet!')
    case 7% FREE DRAINAGE
        dF_dh(P.nz,P.nz) = dz(P.nz)/P.dt.*Ci(P.nz)*P.FrArMtrx(P.nz) + ...
                         + dz(P.nz).*dSINKmacr_dh(P.nz) + ...
                         + P.Kim2(P.nz)./disnod(P.nz) + ...
                         + P.Kip2(P.nz)./(0.5*dz(P.nz)) + ...
                     k*( - dKim2_dhi(P.nz).*hgrad(P.nz) + ...
                         + dKip2_dhi(P.nz) ...
                       );
    otherwise
        error('Wrong setting of W.SwBotB!')
end
%% *ABOVE diagonal
%    *[dF1_dh2]:
dF_dh(1,2)  = - P.Kip2(1)/disnod(2) + ...
              + k*dKip2_dhip1(1)*hgrad(2);
%    [dFi_dhip1(2:P.nz-1,3:P.nz-0)]:
IDX         = sub2ind([P.nz,P.nz],2:P.nz-1,3:P.nz-0);
dF_dh(IDX)  = - P.Kip2(2:P.nz-1)./disnod(3:P.nz-0) + ...
              + k*dKip2_dhip1(2:P.nz-1).*hgrad(3:P.nz-0);
%% *BELOW diagonal
%    [dFi_dhim1(2:P.nz-1,1:P.nz-2)]:
IDX         = sub2ind([P.nz,P.nz],2:P.nz-1,1:P.nz-2);
dF_dh(IDX)  = - P.Kim2(2:P.nz-1)./disnod(2:P.nz-1) + ...
              - k*dKim2_dhim1(2:P.nz-1).*hgrad(2:P.nz-1);
%    *[ dFnz_dhnz-1 ]:
dF_dh(P.nz,P.nz-1) = - P.Kim2(P.nz)/disnod(P.nz) + ...
                     - k*dKim2_dhim1(P.nz)*hgrad(P.nz);
%% clean
clear dz IDX %hgrad disnod dKim2_dhi dKip2_dhi dKim2_dhim1 dKip2_dhip1