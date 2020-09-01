function [Linv, sumTerm] = stehfestAlgorithm(F, t, n)
sumTerm = 0;
for i = 1 : 1 : n
    sumTerm = sumTerm + cSteh(n, i) .* F(i*log(2)/t);
end
Linv = log(2)/t * sumTerm;
end

function C = cSteh(n, i)
sumTerm = 0;
for k = floor( (i+1)/2 ): 1 : min(i, n/2)
   sumTerm = sumTerm + k^(n/2 )*factorial(2*k) / ( factorial(n/2-k) * factorial(k) * ...
       factorial(k-1)*factorial(i-k)*factorial(2*k-i) ); 
end
C = (-1)^(i+n/2)*sumTerm;
end