function [h2]=fnsyst_OLD(j,dt,dz,dzbot,nz,hbot,hsurf,qsurf,qbot,itbc,ibbc,h1,teta,kond,cap,sink,tetas,tetar,alfrs,fi,alfvg,en,alfvg2,en2,ifr,k0,k0macr,bita,bita2,ifc)

%definire km(i) [� il (ki-1/2)] e il kp(i) [� il (ki+1/2)] 

        kp(1)=(kond(1,j)+kond(2,j))/2;
        teta_hbot=fnteta_OLD(hbot,tetas(nz),tetar(nz),alfrs(nz),fi(nz),alfvg(nz),en(nz),alfvg2(nz),en2(nz),ifr(nz));
        teta_hsurf=fnteta_OLD(hsurf,tetas(1),tetar(1),alfrs(1),fi(1),alfvg(1),en(1),alfvg2(1),en2(1),ifr(1));
    
    if or(ifc==1,ifc==3)
        kp(nz)=(kond(nz,j)+fncond_OLD(teta_hbot,tetas(nz),tetar(nz),alfrs(nz),fi(nz),alfvg(nz),en(nz),alfvg2(nz),en2(nz),k0(nz),k0macr(nz),bita(nz),bita2(nz),ifc(nz)))/2;
        km(1)=(kond(1,j)+fncond_OLD(teta_hsurf,tetas(1),tetar(1),alfrs(1),fi(1),alfvg(1),en(1),alfvg2(1),en2(1),k0(1),k0macr(1),bita(1),bita2(1),ifc(1)))/2;
    else 
        kp(nz)=(kond(nz,j)+fncond_OLD(hbot,tetas(nz),tetar(nz),alfrs(nz),fi(nz),alfvg(nz),en(nz),alfvg2(nz),en2(nz),k0(nz),k0macr(nz),bita(nz),bita2(nz),ifc(nz)))/2;
        km(1)=(kond(1,j)+fncond_OLD(hsurf,tetas(1),tetar(1),alfrs(1),fi(1),alfvg(1),en(1),alfvg2(1),en2(1),k0(1),k0macr(1),bita(1),bita2(1),ifc(1)))/2;
    end
        
        km(nz)=(kond(nz,j)+kond(nz-1,j))/2;
    for i=2:(nz-1)
        kp(i)=(kond(i,j)+kond(i+1,j))/2;
        km(i)=(kond(i,j)+kond(i-1,j))/2;
    end

%risolvere il sistema (vedi swap)
%
%intermediate nodes
    for i=2:(nz-1) 
        ratio(i)=dt/(dz(i)^2);
        alfa(i)=-ratio(i)*km(i);
        beta(i)=cap(i,j)+ratio(i)*km(i)+ratio(i)*kp(i);
        gamma(i)=-ratio(i)*kp(i);
        delta(i)=cap(i,j)*h1(i,j)+dz(i)*ratio(i)*(km(i)-kp(i))-dt*sink(i,j);
    end
%
%definire itbc (top boundary condition) e ibbc(bottom boundary condition) (=0 per q e 1 per h) 
%top node flux (positive upward)
    ratio(1)=2*dt/(dz(i)^2);
    if itbc==0     
        alfa(1)=0;
        beta(1)=cap(1,j)+ratio(1)*kp(1);
        gamma(1)=-ratio(1)*kp(1);
        delta(1)=cap(1,j)*h1(1,j)-dz(i)*ratio(1)*(qsurf+kp(1))-dt*sink(1,j);
    else
        alfa(1)=0;
        beta(1)=cap(1,j)+ratio(1)*km(1)+ratio(1)*kp(1);
        gamma(1)=-ratio(1)*kp(1);
        delta(1)=cap(1,j)*h1(1,j)-dz(i)*ratio(1)*(km(1)-kp(1))+ratio(1)*km(1)*hsurf-dt*sink(1,j);
    end
%bottom node flux
    ratio(nz)=2*dt/(dzbot*dz(nz));
    if ibbc==0      
       alfa(nz)=-ratio(nz)*km(nz);
       beta(nz)=cap(nz,j)+ratio(nz)*km(nz);
       gamma(nz)=0;
       delta(nz)=cap(nz,j)*h1(nz,j)+dz(nz)*ratio(nz)*(km(nz)+qbot)-dt*sink(nz,j);
    else
       alfa(nz)=-ratio(nz)*km(nz);
       beta(nz)=cap(nz,j)+ratio(nz)*km(nz)+ratio(nz)*kp(nz);
       gamma(nz)=0;
       delta(nz)=cap(nz,j)*h1(nz,j)+dz(nz)*ratio(nz)*(km(nz)-kp(nz))+ratio(nz)*kp(nz)*hbot-dt*sink(nz,j);
    end

       BB(1)=gamma(1)/beta(1);
       AA(1)=delta(1)/beta(1);
    for i=2:nz
        CC(i)=beta(i)-alfa(i)*BB(i-1);
        AA(i)=(delta(i)-alfa(i)*AA(i-1))/CC(i);
        BB(i)=gamma(i)/CC(i);
    end
        h2(nz)=AA(nz);
    for i=(nz-1):-1:1 
        h2(i)=AA(i)-BB(i)*h2(i+1);
end