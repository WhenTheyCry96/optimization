% To plot selected generation

load('Sinc_3_rand_F.mat');
imprime(1,vxp,vyp,vzp,gen(1,:,1),gen(1,:,2),gen(1,:,3))

figure(2)
plot(optx)
title('Variable Optimization')
xlabel('Iterations')
ylabel('Best Values')
grid on
hold on
plot(opty)
hold off
legend('x_1','x_2')

figure(3)
semilogy(opta)
title('Minimization')
xlabel('Iterations')
ylabel('Best f(x_1,x_2)')
grid on

txt2 = ['f Best: ', num2str(fbest)];
text(0,1,txt2,'Units','normalized',...
     'HorizontalAlignment','left','VerticalAlignment','bottom');

txt3 = ['f Found: ', num2str(optval)];
text(1,1,txt3,'Units','normalized',...
     'HorizontalAlignment','right','VerticalAlignment','bottom');

function imprime(PRINT,vx,vy,vz,x,y,fx)

    if PRINT == 1
        meshc(vx,vy,vz)
        hold on
        title('Maximization')
        xlabel('x')
        ylabel('y')
        zlabel('f(x,y)')
        plot3(x,y,fx,'k*')
        colormap jet
        drawnow
        hold off
        pause(0.1)
    end

end