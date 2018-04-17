function x = heaviside(x)

    idx = x > 0;
    idx1 = x == 0;
    x(idx) = 1;
    x(~idx) = 0;
    x(idx1) = .5;
    