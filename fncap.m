function cap = fncap( x, P, ii )
%     --- calcola la capacita' capillare ---
%     argomenti:
%     x  = potenziale matriciale

df              = P.tetas(ii)-P.tetar(ii);
if x>=0 
   cap          = 0;
   
elseif P.ifr(ii)==1 
    em          = 1-1/P.en(ii);
    cap         = df*(1+abs(P.alfvg(ii)*x)^P.en(ii))^(-em-1)*em * ...
                  abs(P.alfvg(ii)*x)^(P.en(ii)-1)*P.alfvg(ii)*P.en(ii);
    
elseif P.ifr(ii)==2
    em2         = 1-1/P.en2(ii);
    piece1      = -P.alfrs(ii)^2*abs(x)*exp(-P.alfrs(ii)*abs(x));
    piece2      = P.alfvg2(ii)*P.en2(ii)*em2 * ...
                  (P.alfvg2(ii)*abs(x))^(P.en2(ii)-1) * ...
                  (1+(P.alfvg2(ii)*abs(x))^P.en2(ii))^(-1-em2);
    cap         = df*(P.fi(ii)*piece1 + (1-p.fi(ii))*piece2);

elseif p.ifr(ii)==3
    em          = 1-1/P.en(ii);
    em2         = 1-1/P.en2(ii);
    piece1      = P.alfvg(ii)*P.en(ii)*em * ...
                  (P.alfvg(ii)*abs(x))^(P.en(ii)-1) * ...
                  (1+(P.alfvg(ii)*abs(x))^P.en(ii))^(-1-em);
    piece2      = P.alfvg2(ii)*P.en2(ii)*em2 * ...
                  (P.alfvg2(ii)*abs(x))^(P.en2(ii)-1) * ...
                  (1+(P.alfvg2(ii)*abs(x))^P.en2(ii))^(-1-em2);
    cap         = df*(P.fi(ii)*piece1 + (1-P.fi(ii))*piece2);
end

return