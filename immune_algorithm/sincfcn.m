function z = sincfcn(xx)

    % Sinc function

    % Search domain: [0,20]
    % Global maximum: f(xx) = 1 | x = (10,10)

    % xx = max(0,min(20,xx));

    z = sin(sum(abs(xx-10),2))./sum(abs(xx-10),2);

end