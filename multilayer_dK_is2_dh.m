function dK_is2_dh = multilayer_dK_is2_dh( dK_dh, method, K_is1, K_i, dz_is1, dz_i, im1 )
% dK_is2_dh = multilayer_dK_is2_dh( dK_dh, method, K_is1, K_i, dz_is1, dz_i, im1 )
% 
% NOTE
%   This function is very similar to the "dkmean" function of SWAP used in
%   headcalc.for and implemented in headcalc.for.
% 
% DESCRIPTION
%   This function computes the derivative of the internodal conductivity to
%   the pressure head at current i node using an implicit linearization of
%   hydraulic conductivity.
% 
% INPUTs
%   dK_dh:      Derivative of the conductivity to the pressure head at
%               specific node. It can assume one of the following forms:
%                  -----------------------------------------------
%                   nodal der. conduc. | internodal der. conduc.
%                        [dK_dh]       |       [dK_is2_dh]
%                  -----------------------------------------------
%                   *dK(i-1) / dh(i-1) | if dK(i-1/2) / dh(i-1)
%                   *dK(i)   / dh(i)   | if dK(i-1/2) / dh(i)
%                   *dK(i)   / dh(i)   | if dK(i+1/2) / dh(i)
%                   *dK(i+1) / dh(i+1) | if dK(i+1/2) / dh(i+1)
%                  -----------------------------------------------
% 
%   method:     One of the following methods:
%                   *arithmic mean
%                   *weighted arithmic mean
%                   *geometric mean
%                   *weighted geometric mean
% 
%   K_is1:      It is the hydraulic conductivity at one of the following
%               nodes ("s" stays for "+" or "-"):
%                   *(i-1)  and "s" stays for i "minus" 1
%                   *(i+1)  and "s" stays for i "plus"  1
% 
%   K_i:        It is the hydraulic conductivity at node i.
% 
%   dz_is1:     It is the compartiment thickness at one of the following
%               nodes ("s" stays for "+" or "-"):
%                   *(i-1)  and "s" stays for i "minus" 1
%                   *(i+1)  and "s" stays for i "plus"  1
% 
%   dz_i:       It is the compartiment thickness at node i.
% 
%   im1:        Flag to set whether the derivative of the internodal
%               conductivity (at i-1/2 or at i+1/2) must be calculated to
%               the pressure head at node i-1 (im1=true) or at node i
%               (im1=false).
%                   *true   --> dK_dh was computed as dK(i-1)/dh(i-1) or
%                               dK(i+1)/dh(i+1)
%                   *false  --> dK_dh was computed as dK(i)  /dh(i)
%               SWAP-32 instead of a flag such as im1 calls the same
%               function inverting the order of K's and dz's (compare for
%               instance the contrasting computation of dkmean at lines 333
%               and 336!!).
% 
% OUTPUTs
%   dK_is2_dh:  Derivative of the internodal conductivity to the pressure
%               head.
%% main
switch method
    case 1% arithmic mean
        % im1 does not affect this case!
        dK_is2_dh = dK_dh/2;
    case 2% weighted arithmic mean
        if im1
            dK_is2_dh = dz_is1  ./ (dz_is1 + dz_i) .* dK_dh;
        else
            dK_is2_dh = dz_i    ./ (dz_is1 + dz_i) .* dK_dh;
        end
    case 3% geometric mean
        if im1
            dK_is2_dh = (K_i./K_is1) .* dK_dh /2;
        else
            dK_is2_dh = (K_is1./K_i) .* dK_dh /2;
        end
    case 4% weighted geometric mean
        if im1
            dK_is2_dh = dz_is1./(dz_is1 + dz_i) .* (K_i./K_is1).^(dz_i./(dz_is1+dz_i))  .* dK_dh;
        else
            dK_is2_dh = dz_i./(dz_is1 + dz_i)   .* (K_is1./K_i).^(dz_is1./(dz_is1+dz_i)).* dK_dh;
        end
    otherwise
        error('You defined a wrong method to calculate the internodal conductivity')
end
%% end
return