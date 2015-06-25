function r = hypergeometricRatio(D, x)
   r = hypergeom(2, D+1, x) ./ (hypergeom(1, D, x) * D);
end

