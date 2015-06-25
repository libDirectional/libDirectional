function r = gradHypergeometricRatio(D, x)
    r = 2/(D+1)/D * hypergeom(3, D+2, x) ./ hypergeom(1, D, x) ...
        - 1/D^2 * hypergeom(2, D+1, x).^2 ./ hypergeom(1, D, x).^2;
end

