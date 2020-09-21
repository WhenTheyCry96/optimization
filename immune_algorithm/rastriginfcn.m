function z = rastriginfcn(xx)

    % Rastrigin's function

    % Search domain: [-5.12,5.12]
    % Global minimum: f(x) = 0 | x = (0,...,0)

    % xx = max(-5.12,min(5.12,xx));

    z = 10.0*size(xx,2) + sum(xx.^2 - 10.0*cos(2*pi.*xx),2);

end