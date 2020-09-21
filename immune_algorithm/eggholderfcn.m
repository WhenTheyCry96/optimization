function z = eggholderfcn(xx)

    % Eggholder function

    % Search domain: [-512,512]
    % Global minimum: f(x) = -959.6407 | x = (512,404.2319)

    % xx = max(-512,min(512,xx));

    x1 = xx(:,1);
    x2 = xx(:,2);

    term1 = - (x2 + 47).*sin(sqrt(abs(x2 + x1/2 + 47)));
    term2 = - x1.*sin(sqrt(abs(x1 - (x2 + 47))));

    z = -term1 - term2;

end