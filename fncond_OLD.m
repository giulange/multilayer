function [kond]= fncond(x,tetas,tetar,alfrs,fi,alfvg,en,alfvg2,en2,k0,k0macr,bita,bita2,ifc)
% x=teta=contenuto d'acqua (ifc=1;3)
% x=h=potenziale (ifc=2 (DURNER); ifc=4 (R&S interacting distributions); ifc=5 (R&S independent distributions [see R&S 1993 eq.21])

% % global x_global
% % % global bimwr
% % global lx 
% % global ly


if ifc==1
    em=1-1/en;
    sef=(x-tetar)/(tetas-tetar);
    if sef>0.99999   
        kond=k0;
    else 
        kond=k0*(sef)^bita*(1.-(1.-sef^(1/em))^em)^2;
    end 
    
elseif ifc==2
    em=1-1/en;
    em2=1-1/en2;    
    sef1=(1+(alfvg*abs(x))^en)^-em;
    sef2=(1+(alfvg2*abs(x))^en2)^-em2;
    sef=fi*sef1+(1-fi)*sef2;
    if sef>0.99999   
        kond=k0;
    else 
        kr_macr=fi*alfvg*(1-(1-sef1^(1/em))^em);
        kr_micr=(1-fi)*alfvg2*(1-(1-sef2^(1/em2))^em2);
        kond=k0*(fi*sef1+(1-fi)*sef2)^bita*((kr_macr+kr_micr)/(fi*alfvg+(1-fi)*alfvg2))^2;
    end
elseif ifc==3
    sef=(x-tetar)/(tetas-tetar);
    if sef>0.99999
        kond=k0;
    else 
        kond=k0*(fi*exp(-bita*(tetas-x))+(1.-fi)*exp(-bita2*(tetas-x)));
    end
    
elseif ifc==4
    em2=1-1/en2;
    sef1=(1+alfrs*abs(x))*exp(-alfrs*abs(x));
    sef2=(1+(alfvg2*abs(x))^en2)^-em2;
    sef=fi*sef1+(1-fi)*sef2;
    if sef>0.99999   
        kond=k0;
    else 
        gmacr=alfrs*exp(-alfrs*abs(x));        
        gmicr=alfvg2*en2*betainc(sef2^(1/em2),1,em2);        
        gmacr_0=alfrs;
        gmicr_0=alfvg2*en2;
        kond=k0*sef^bita*(fi*gmacr+(1-fi)*gmicr)/(fi*gmacr_0+(1-fi)*gmicr_0);
    end
    
elseif ifc==5
    em2=1-1/en2;
    sef1=(1+alfrs*abs(x))*exp(-alfrs*abs(x));
    sef2=(1+(alfvg2*abs(x))^en2)^-em2;
    sef=fi*sef1+(1-fi)*sef2;
    if sef>0.99999   
        kond=k0macr;
    else 
        kr_macr=sef1^bita*exp(-2*alfrs*abs(x));
        kr_micr=sef2^bita*(1.-(1.-sef2^(1/em2))^em2)^2;
        kond=k0macr*kr_macr+k0*kr_micr;
    end    
end