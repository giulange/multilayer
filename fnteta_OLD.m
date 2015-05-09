function [teta]=fnteta(x,tetas,tetar,alfrs,fi,alfvg,en,alfvg2,en2,ifr)
%x=h

df=tetas-tetar;

if x>0 
    teta=tetas;
elseif ifr==1
    em=1-1/en;
    teta=tetar+df/(1.+(abs(x)*alfvg)^en)^em;

elseif ifr==2
    em2=1-1/en2;
    teta=tetar+df*(fi*((1+alfrs*abs(x))*exp(-alfrs*abs(x)))+(1-fi)*(1+(alfvg2*abs(x))^en2)^(-em2));
    
%     lx=-1000000:1000:0;
%     ly=tetar+df.*(fi.*((1+alfrs.*abs(lx)).*exp(-alfrs.*abs(lx)))+(1-fi).*(1+(alfvg.*abs(lx)).^en).^(-em));
    
elseif ifr==3
    em=1-1/en;
    em2=1-1/en2;    
    sef1=(1+(alfvg*abs(x))^en)^-em;
        sef2=(1+(alfvg2*abs(x))^en2)^-em2;
    sef=fi*sef1+(1-fi)*sef2;
    teta=tetar+df*sef;    
end

