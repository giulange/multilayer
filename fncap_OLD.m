function [cap]=fncap(x,tetas,tetar,alfrs,fi,alfvg,en,alfvg2,en2,ifr)
%     --- calcola la capacita' capillare ---
%     argomenti:
%     x  = potenziale matriciale

df=tetas-tetar;
if x>=0 
   cap=0;
   
elseif ifr==1 
    em=1-1/en;
    cap=df*(1+abs(alfvg*x)^en)^(-em-1)*em*abs(alfvg*x)^(en-1)*alfvg*en;
    
elseif ifr==2
    em2=1-1/en2;
    piece1=-alfrs^2*abs(x)*exp(-alfrs*abs(x));  
    piece2=alfvg2*en2*em2*(alfvg2*abs(x))^(en2-1)*(1+(alfvg2*abs(x))^en2)^(-1-em2);
    cap=df*(fi*piece1+(1-fi)*piece2);
elseif ifr==3
    em=1-1/en;
    em2=1-1/en2;
    piece1=alfvg*en*em*(alfvg*abs(x))^(en-1)*(1+(alfvg*abs(x))^en)^(-1-em);
    piece2=alfvg2*en2*em2*(alfvg2*abs(x))^(en2-1)*(1+(alfvg2*abs(x))^en2)^(-1-em2);
    cap=df*(fi*piece1+(1-fi)*piece2);
end
