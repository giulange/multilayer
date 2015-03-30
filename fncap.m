function cap = fncap( x, Psh, ii )
% cap = fncap( x, P.sh, ii )
%     --- calcola la capacita' capillare ---
%     argomenti:
%     x  = potenziale matriciale

df              = Psh.tetas(ii)-Psh.tetar(ii);
if x>=0 
   cap          = 0;
   
elseif Psh.ifr(ii)==1 
    em          = 1-1/Psh.en(ii);
    cap         = df*(1+abs(Psh.alfvg(ii)*x)^Psh.en(ii))^(-em-1)*em * ...
                  abs(Psh.alfvg(ii)*x)^(Psh.en(ii)-1)*Psh.alfvg(ii)*Psh.en(ii);
    
elseif Psh.ifr(ii)==2
    em2         = 1-1/Psh.en2(ii);
    piece1      = -Psh.alfrs(ii)^2*abs(x)*exp(-Psh.alfrs(ii)*abs(x));
    piece2      = Psh.alfvg2(ii)*Psh.en2(ii)*em2 * ...
                  (Psh.alfvg2(ii)*abs(x))^(Psh.en2(ii)-1) * ...
                  (1+(Psh.alfvg2(ii)*abs(x))^Psh.en2(ii))^(-1-em2);
    cap         = df*(Psh.fi(ii)*piece1 + (1-Psh.fi(ii))*piece2);

elseif p.ifr(ii)==3
    em          = 1-1/Psh.en(ii);
    em2         = 1-1/Psh.en2(ii);
    piece1      = Psh.alfvg(ii)*Psh.en(ii)*em * ...
                  (Psh.alfvg(ii)*abs(x))^(Psh.en(ii)-1) * ...
                  (1+(Psh.alfvg(ii)*abs(x))^Psh.en(ii))^(-1-em);
    piece2      = Psh.alfvg2(ii)*Psh.en2(ii)*em2 * ...
                  (Psh.alfvg2(ii)*abs(x))^(Psh.en2(ii)-1) * ...
                  (1+(Psh.alfvg2(ii)*abs(x))^Psh.en2(ii))^(-1-em2);
    cap         = df*(Psh.fi(ii)*piece1 + (1-Psh.fi(ii))*piece2);
end

return