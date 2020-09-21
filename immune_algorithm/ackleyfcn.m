function z = ackleyfcn(xx)

    % Ackley's function

    % Search domain: [-32,32]
    % Global minimum: f(x) = 0 | x = (0,...,0)

    % xx = max(-32,min(32,xx));

    d = size(xx,2);

    z = -20*exp(-0.2*sqrt(1/d*sum(xx.^2,2))) - exp(1/d*sum(cos(2*pi*xx),2)) + 20 + exp(1);

end