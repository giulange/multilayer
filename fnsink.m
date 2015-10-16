function sink = fnsink( x, z_i,op,ptra,Droot,ifg,ifs,V)
% fnsink( x, dpt,op,ptra,Droot,ifg,ifs,V)
% 
% DESCRIPTION
%   This function computes root water extraction flux(??) for each soil
%   compartment. It computes the sink for the current compartment passed as
%   input.
% 
% INPUTs
%   x:          Pressure head of soil compartment, e.g. h_jm1(P.nz)
%   dpt:        e.g. z(i)
%   op:         e.g. EC(i,j)
% 
% OUTPUTs
%   sink:       ??
% 
% NOTEs
%   These are unused but they were passed to old function: V.hIIIH

%% null sink (Droot is null)
if Droot<=0
    sink            = 0;
    return
end
%% PRE
arw = 0;
ars = 0;
%% ifg (??)
switch ifg% {1,2,3,4}
    % distribuzione radici
    case 1% linear distribution (constant with depth)
        gz          = 1/Droot;
        
    case 2% logistic distribution
        gz          = V.rda*V.rdb*V.rdc*exp(V.rdc*z_i) / ...
                      (V.rdb+exp(V.rdc*z_i))^2;
                  
    case 3% Prasad-type distribution (triagle)
        gz          = 2*(Droot-z_i)/Droot^2;

    case 4% two-linear distribution
        if Droot<=V.zc
            gz      = 1/Droot;
        else
           A1f      = V.g0*V.zc;
           A2f      = V.gzc*(V.Drf-V.zc);
           A3t      = V.gzc*(V.Drf-Droot);
           delt_A1  = A1f*A3t;
           delt_A2  = A2f*A3t;
           A1t      = A1f+delt_A1;
           A2t      = A2f-A3t+delt_A2;
           if z_i<V.zc
               gz   = A1t/V.zc;
           else
               gz   = A2t/(Droot-V.zc);
           end
        end

end
%% ifs (??)
Sa                  = ptra*gz; % sink potenziale
switch ifs% quale sink vuoi usare?? --> {1,2,3,4,5,6,7}
    
    % FEDDES water reduction factor
    case 1
        if x>=V.hI
            arw     = 0;% fattore di riduzione
        elseif and(x>=V.hII,x<V.hI)
            arw     = (V.hI-x)/(V.hI-V.hII);
        elseif and(x>=V.hIIIL,x<V.hII)
            arw     = 1;
        elseif and(x>=V.hIV,x<V.hIIIL)
            arw     = (x-V.hIV)/(V.hIIIL-V.hIV);
        end
        sink        = Sa*arw;% sink effettivo che entra in Richards

    % FEDDES water reduction factor
    case 2
        if x>=V.hI
            arw     = 0;
        elseif and(x>=V.hII,x<V.hI)
            arw     = (V.hI-x)/(V.hI-V.hII);
        elseif and(x>=V.hIIIL,x<V.hII)
            arw     = 1;
        elseif and(x>=V.hIV,x<V.hIIIL)
            arw     = (x-V.hIV)/(V.hIIIL-V.hIV);
        end
        sink        = Sa*arw;
    
    % van Genuchten water reduction factor
    case 3
        arw         = 1/(1+(x/V.hw50)^V.pw1);
        sink        = Sa*arw;
    
    % Maas&Hoffman salinity reduction factor
    case 4
        if and(op>=V.aMH,op<=0)
            ars     = 1;
        elseif and(op>=V.aMH-1/V.bMH,op<=V.aMH)
            ars     = 1+V.bMH*(op-V.aMH);
        elseif op<=V.aMH-1/V.bMH
            ars     = 0;
        end
        sink        = Sa*ars;

    % van Genuchten salinity reduction factor
    % distribution
    case 5
        Sa          = ptra*gz;
        ars         = 1/(1+(op/V.hs50)^V.ps1);
        sink        = Sa*ars;

    % van Genuchten multiplicative water and salinity reduction factor with
    % uniform root distribution
    case 6
        Sa          = ptra*gz;
        arw         = 1/(1+(x/V.hw50)^V.pw1);
        ars         = 1/(1+(op/V.hs50)^V.ps1);
        sink        = Sa*arw*ars;

    % Multiplicative FEDDES water and Maas&Hoffman salinity reduction
    % factors with uniform root distribution
    case 7
        Sa          = ptra*gz;
        if x>=V.hI
            arw     = 0;
        elseif and(x>=V.hII,x<V.hI)
            arw     = (V.hI-x)/(V.hI-V.hII);
        elseif and(x>=V.hIIIL,x<V.hII)
            arw     = 1;
        elseif and(x>=V.hIV,x<V.hIIIL)
            arw     = (x-V.hIV)/(V.hIIIL-V.hIV);
        end
        if and(op>=V.aMH,op<=0)
            ars     = 1;
        elseif and(op>=V.aMH-1/V.bMH,op<=V.aMH)
            ars     = 1+V.bMH*(op-V.aMH);
        elseif op<=V.aMH-1/V.bMH
            ars     = 0;
        end
        sink        = Sa*arw*ars;
end
%% end
return