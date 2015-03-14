function sink = fnsink( x, P, W, V )
%x=h1(i,j)
%dpt=z(i)
%op=EC(i,j)

% unused but passed to old function:
%   V.hIIIH, W.nz

arw                     = 0;
ars                     = 0;
if P.Droot>0
    switch V.ifg
        case 1
            gz          = 1/P.Droot;
        case 2
            gz          = V.rda*V.rdb*V.rdc*exp(V.rdc*P.dpt) / ...
                          (V.rdb+exp(V.rdc*P.dpt))^2;
        case 3
            gz          = 2*(P.Droot-P.dpt)/P.Droot^2;
            
        case 4
            if P.Droot<=V.zc
                gz      = 1/P.Droot;
            else
               A1f      = V.g0*V.zc;
               A2f      = V.gzc*(V.Drf-V.zc);
               A3t      = V.gzc*(V.Drf-P.Droot);
               delt_A1  = A1f*A3t;
               delt_A2  = A2f*A3t;
               A1t      = A1f+delt_A1;
               A2t      = A2f-A3t+delt_A2;
               if P.dpt<V.zc
                   gz   = A1t/V.zc;
               else
                   gz   = A2t/(P.Droot-V.zc);
               end
            end
    end
end

% FEDDES water reduction factor with uniform root distribution
if V.ifs==1
    if P.Droot>0
        Sa              = W.tp*gz;
        if x>=V.hI
            arw         = 0;
        elseif and(x>=V.hII,x<V.hI)
            arw         = (V.hI-x)/(V.hI-V.hII);
        elseif and(x>=V.hIIIL,x<V.hII)
            arw         = 1;
        elseif and(x>=V.hIV,x<V.hIIIL)
            arw         = (x-V.hIV)/(V.hIIIL-V.hIV);
        end
    sink                = Sa*arw;
    else
    sink=0;
    end
end

% FEDDES water reduction factor with logistic root distribution
if V.ifs==2
    if P.Droot>0
    Sa                  = W.tp*gz;
        if x>=V.hI
            arw         = 0;
        elseif and(x>=V.hII,x<V.hI)
            arw         = (V.hI-x)/(V.hI-V.hII);
        elseif and(x>=V.hIIIL,x<V.hII)
            arw         = 1;
        elseif and(x>=V.hIV,x<V.hIIIL)
            arw         = (x-V.hIV)/(V.hIIIL-V.hIV);
        end
        sink            = Sa*arw;
    
    else
        sink            = 0;
    end
end
    
% van Genuchten water reduction factor with uniform root distribution
if V.ifs==3
    if P.Droot>0
        Sa              = W.tp*gz;
        arw             = 1/(1+(x/V.hw50)^V.pw1);
        sink            = Sa*arw;
    else
        sink            = 0;
    end
end
    
% Maas&Hoffman salinity reduction factor with uniform root distribution
if V.ifs==4
    if P.Droot>0
    Sa                  = W.tp*gz;
        if and(P.op>=V.aMH,P.op<=0)
            ars         = 1;
        elseif and(P.op>=V.aMH-1/V.bMH,P.op<=V.aMH)
            ars         = 1+V.bMH*(P.op-V.aMH);
        elseif P.op<=V.aMH-1/V.bMH
            ars         = 0;
        end
        sink            = Sa*ars;
    else
        sink            = 0;
    end
end

% van Genuchten salinity reduction factor with uniform root distribution
if V.ifs==5
    if P.Droot>0
        Sa              = W.tp*gz;
        ars             = 1/(1+(P.op/V.hs50)^V.ps1);
        sink            = Sa*ars;
    else
        sink            = 0;
    end
end

% van Genuchten multiplicative water and salinity reduction factor with uniform root distribution
if V.ifs==6
    if P.Droot>0
        Sa              = W.tp*gz;
        arw             = 1/(1+(x/V.hw50)^V.pw1);
        ars             = 1/(1+(P.op/V.hs50)^V.ps1);
        sink            = Sa*arw*ars;
    else
        sink            = 0;
    end
end

% Multiplicative FEDDES water and Maas&Hoffman salinity reduction factors with uniform root distribution
if V.ifs==7
    if P.Droot>0
        Sa              = W.tp*gz;
        if x>=V.hI
            arw         = 0;
        elseif and(x>=V.hII,x<V.hI)
            arw         = (V.hI-x)/(V.hI-V.hII);
        elseif and(x>=V.hIIIL,x<V.hII)
            arw         = 1;
        elseif and(x>=V.hIV,x<V.hIIIL)
            arw         = (x-V.hIV)/(V.hIIIL-V.hIV);
        end
        if and(P.op>=V.aMH,P.op<=0)
            ars         = 1;
        elseif and(P.op>=V.aMH-1/V.bMH,P.op<=V.aMH)
            ars         = 1+V.bMH*(P.op-V.aMH);
        elseif P.op<=V.aMH-1/V.bMH
           ars          = 0;
        end
        sink            = Sa*arw*ars;
    else
        sink            = 0;
    end
end

return