function z = bealefcn(xx)

    % Beale's function

    % Search domain: [-4.5,4.5]
    % Global minimum: f(xx) = 0 | xx = (3,0.5)

    % xx = max(-5.12,min(5.12,xx));
    
    x1 = xx(:,1);
    x2 = xx(:,2);

    z = (1.5 - x1 + x1.*x2).^2 + (2.25 - x1 + x1.*x2.^2).^2 + (2.625 - x1 + x1.*x2.^3).^2;

end