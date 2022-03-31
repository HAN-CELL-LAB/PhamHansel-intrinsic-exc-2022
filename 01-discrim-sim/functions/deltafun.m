function Y = deltafun(X, epsi_deltafun)
    Y = 1.0 * (abs(X) < epsi_deltafun);
end
