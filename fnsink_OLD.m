function [sink]=fnsink(x,dpt,op,Tp,hI,hII,hIIIH,hIIIL,hIV,hw50,pw1,hs50,ps1,aMH,bMH,Droot,ifs,rda,rdb,rdc,nz,ifg,zc,g0,gzc,Drf)
%x=h1(i,j)
%dpt=z(i)
%op=EC(i,j)

% input('sono in sink')

arw=0;
ars=0;
if Droot>0
    if ifg==1
        gz=1/Droot;
    elseif ifg==2
        gz=rda*rdb*rdc*exp(rdc*dpt)/(rdb+exp(rdc*dpt))^2;
    elseif ifg==3
        gz=2*(Droot-dpt)/Droot^2;
    elseif ifg==4
        if Droot<=zc
            gz=1/Droot;
        else
           A1f=g0*zc;
           A2f=gzc*(Drf-zc);
           A3t=gzc*(Drf-Droot);
           delt_A1=A1f*A3t;
           delt_A2=A2f*A3t;
           A1t=A1f+delt_A1;
           A2t=A2f-A3t+delt_A2;
           if dpt<zc
               gz=A1t/zc;
           else
               gz=A2t/(Droot-zc);
           end
        end
    end
end



% FEDDES water reduction factor with uniform root distribution
    if ifs==1
        if Droot>0
        Sa=Tp*gz;
            if x>=hI
            arw=0;
            elseif and(x>=hII,x<hI)
            arw=(hI-x)/(hI-hII);
            elseif and(x>=hIIIL,x<hII)
            arw=1;
            elseif and(x>=hIV,x<hIIIL)
            arw=(x-hIV)/(hIIIL-hIV);
            end
        sink=Sa*arw;
        else
        sink=0;
        end
    end
    
 % FEDDES water reduction factor with logistic root distribution
    if ifs==2
        if Droot>0
        Sa=Tp*gz;
            if x>=hI
            arw=0;
            elseif and(x>=hII,x<hI)
            arw=(hI-x)/(hI-hII);
            elseif and(x>=hIIIL,x<hII)
            arw=1;
            elseif and(x>=hIV,x<hIIIL)
            arw=(x-hIV)/(hIIIL-hIV);
            end
        sink=Sa*arw;
        else
        sink=0;
        end
    end
    
    % van Genuchten water reduction factor with uniform root distribution
    if ifs==3
        if Droot>0
        Sa=Tp*gz;
        arw=1/(1+(x/hw50)^pw1);
        sink=Sa*arw;
        else
        sink=0;
        end
    end
    
    % Maas&Hoffman salinity reduction factor with uniform root distribution
    if ifs==4
        if Droot>0
        Sa=Tp*gz;
            if and(op>=aMH,op<=0)
            ars=1;
            elseif and(op>=aMH-1/bMH,op<=aMH)
            ars=1+bMH*(op-aMH);
            elseif op<=aMH-1/bMH
            ars=0;
            end
        sink=Sa*ars;
        else
        sink=0;
        end
    end
    
        % van Genuchten salinity reduction factor with uniform root distribution
    if ifs==5
        if Droot>0
        Sa=Tp*gz;
        ars=1/(1+(op/hs50)^ps1);
        sink=Sa*ars;
        else
        sink=0;
        end
    end
    
        % van Genuchten multiplicative water and salinity reduction factor with uniform root distribution
    if ifs==6
        if Droot>0
        Sa=Tp*gz;
        arw=1/(1+(x/hw50)^pw1);
        ars=1/(1+(op/hs50)^ps1);
        sink=Sa*arw*ars;
        else
        sink=0;
        end
    end
    
    % Multiplicative FEDDES water and Maas&Hoffman salinity reduction factors with uniform root distribution
    if ifs==7
        if Droot>0
        Sa=Tp*gz;
            if x>=hI
            arw=0;
            elseif and(x>=hII,x<hI)
            arw=(hI-x)/(hI-hII);
            elseif and(x>=hIIIL,x<hII)
            arw=1;
            elseif and(x>=hIV,x<hIIIL)
            arw=(x-hIV)/(hIIIL-hIV);
            end

            if and(op>=aMH,op<=0)
            ars=1;
            elseif and(op>=aMH-1/bMH,op<=aMH)
            ars=1+bMH*(op-aMH);
            elseif op<=aMH-1/bMH
            ars=0;
            end
        sink=Sa*arw*ars;
        else
        sink=0;
        end
    end



