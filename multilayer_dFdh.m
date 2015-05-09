%% NOTES
% This script computes ...

%% init
dKim2_dhi   = zeros(P.nz,1);
dKip2_dhi   = zeros(P.nz,1);
dKim2_dhim1 = zeros(P.nz,1);
dKip2_dhip1 = zeros(P.nz,1);
% 
dF_dh       = zeros(P.nz,P.nz);
%% *intermediate nodes [diagonal, below-diagonal,above-diagonal :: 2:P.nz-1]
% UPDATE :: implicit linearization
% Derivative of the internodal hydraulic conductivity to the pressure
% head (SWAP-32 manual, Appendix 5, page 239):
dKim2_dhi(i_i)   = multilayer_dK_is2_dh( dKi_dhi(i_i),   W.Kmeth, K(i_im1), K(i_i), dz_im1, dz_i, false );
dKip2_dhi(i_i)   = multilayer_dK_is2_dh( dKi_dhi(i_i),   W.Kmeth, K(i_ip1), K(i_i), dz_ip1, dz_i, false );
dKim2_dhim1(i_i) = multilayer_dK_is2_dh( dKi_dhi(i_im1), W.Kmeth, K(i_im1), K(i_i), dz_im1, dz_i, true  );
dKip2_dhip1(i_i) = multilayer_dK_is2_dh( dKi_dhi(i_ip1), W.Kmeth, K(i_ip1), K(i_i), dz_ip1, dz_i, true  );
% top    node particular cases:     check with Antonio
dKip2_dhi(1)     = multilayer_dK_is2_dh( dKi_dhi(1), W.Kmeth, K(2), K(1), dz(2), dz(1), false );
dKip2_dhip1(1)   = multilayer_dK_is2_dh( dKi_dhi(2), W.Kmeth, K(2), K(1), dz(2), dz(1), true  );
% bottom node particular cases:     check with Antonio
dKim2_dhi(P.nz)  = multilayer_dK_is2_dh( dKi_dhi(P.nz),   W.Kmeth, K(P.nz-1), K(P.nz), dz(P.nz-1), dz(P.nz), false );
dKim2_dhim1(P.nz)= multilayer_dK_is2_dh( dKi_dhi(P.nz-1), W.Kmeth, K(P.nz-1), K(P.nz), dz(P.nz-1), dz(P.nz), true  );

% partial derivative of the macro-pore exchange to the pressure head:
dMacrPor    = zeros(P.nz,1); % must be developed!!!

%% Jacobian coefficients of the tri-diagonal system of equations
% (SWAP-32 manual, Appendix 4):
%   -Partial derivatives of Fi to pressure heads parameterised upon k
%   to allow for treating the hydraulic conductivities implicitly (when
%   k=1)-

%   *top     boundary:
%       -dF/dh (derivative of Fi to the pressure head)
%       -OLD :: W.itbc (1-->(~isflux), 0-->isflux)
if ~isflux%     HEAD controlled top boundary
    dK1_2_dhi = 0;  % DEFINE!! ask Antonio
    K1_2 = 0;       % DEFINE!! ask Antonio
    dF_dh(1,1)  = dz(1)/P.dt.*Ci(1) + dz(1).*dMacrPor(1) + ...
          + K1_2/(dz(1)/2) + 2*Kim2(1)./(dz(1)+dz(2)) + ...
          - dK1_2_dhi*(2*(hsurf-h_pm1(1))/(dz(1)/2) +1) + ...
          + k*dKip2_dhi(1).*(2*(h_pm1(1)-h_pm1(2))./(dz(1)+dz(2)) +1);
elseif isflux%  FLUX controlled top boundary
    dF_dh(1,1)  = dz(1)/P.dt.*Ci(1) + dz(1).*dMacrPor(1) + ...
          + 2*Kim2(1)./(dz(1)+dz(2)) + ...
          + k*dKip2_dhi(1).*(2*(h_pm1(1)-h_pm1(2))./(dz(1)+dz(2)) +1);
end

%   *bottom  boundary:
dF_dh(P.nz,P.nz) = dFdh_botB;

%    *diagonal [dFi_dhi(i)]:
IDX         = sub2ind([P.nz,P.nz],i_i,i_i);
dF_dh(IDX)  = dz_i/P.dt.*Ci(i_i)   + dz_i.*dMacrPor(i_i) + ...
          + 2*Kim2(i_i)./(dz_im1+dz_i) + 2*Kip2(i_i)./(dz_i+dz_ip1) + ...
          - k*dKim2_dhi(i_i).*(2*(h_pm1(i_im1)-h_pm1(i_i))./(dz_im1+dz_i) +1) + ...
          + k*dKip2_dhi(i_i).*(2*(h_pm1(i_i)-h_pm1(i_ip1))./(dz_i+dz_ip1) +1);

%    *below diagonal [dFi_dhim1(i)]:
IDX         = sub2ind([P.nz,P.nz],i_i,i_im1);
dF_dh(IDX)  = - 2*Kim2(i_i)./(dz_im1+dz_i) - ...
          k*dKim2_dhim1(i_i).*(2*(h_pm1(i_im1)-h_pm1(i_i))./(dz_im1+dz_i) +1);

%    *above diagonal [dFi_dhip1(i)]:
IDX         = sub2ind([P.nz,P.nz],i_i,i_ip1);
dF_dh(IDX)  = - 2*Kip2(i_i)./(dz_i+dz_ip1) + ...
          k*dKip2_dhip1(i_i).*(2*(h_pm1(i_i)-h_pm1(i_ip1))./(dz_i+dz_ip1) +1);
      
%    *below diagonal [ dFn_dhn-1 ]:
dF_dh(P.nz,P.nz-1) = - 2*Kim2(P.nz)/(dz(P.nz-1)+dz(P.nz)) - ...
            k*dKim2_dhim1(P.nz)*(2*(h_pm1(P.nz-1)-h_pm1(P.nz))/(dz(P.nz-1)+dz(P.nz)) +1);
        
%    *above diagonal [ dF1_dh2 ]:
dF_dh(1,2)  = - 2*Kip2(1)/(dz(1)+dz(2)) + ...
            k*dKip2_dhip1(1)*(2*(h_pm1(1)-h_pm1(2))/(dz(1)+dz(2)) +1);
%