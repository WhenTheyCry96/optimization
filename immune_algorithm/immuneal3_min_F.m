close all
clear all
clc

format long

tic

rng('default')

%% Function Setting
f = @rastriginfcn; %CHANGE ME TO TRY ANOTHER TEST FUNCTION
fbest = 0; % Global f best (minimum), %CHANGE ME ACCORDING TO TEST FUNCTION 

% side constraints
varMin= -5.12; %CHANGE ME ACCORDING TO TEST FUNCTION
varMax= 5.12; %CHANGE ME ACCORDING TO TEST FUNCTION


numvar=2; %CHANGE ME ACCORDING TO TEST FUNCTION

% hyperparameters
Converge=1e-3; %0.00001; % diff [%] %CHANGE ME ACCORDING TO TEST FUNCTION
crovrate=0;
mutrate=0.5;
surv=0.1;
maxit=1000;
search=0.1;

N=50;
mem=zeros(N,2);
xbaby=ones(N,1);
ybaby=ones(N,1);
c=ones(N,1);
e=ones(N,1);

%make grid
gridN=200;
a = meshgrid(linspace(varMin, varMax, gridN));
b = meshgrid(linspace(varMin, varMax, gridN))';
vxp = a;
vyp = b;
vzp = f([a(:),b(:)]);
vzp = reshape(vzp,size(a));

% DATA=[population size, x, y, f(x,y), # of iterations, # of function calls]
% DATA=zeros(100,6);
% for re=1:100


newmem=zeros(N,3);
optx=[];
opty=[];
opta=[];

% initial pop
for i=1:N
    mem(i,1)=varMin+(varMax-varMin)*rand;
    mem(i,2)=varMin+(varMax-varMin)*rand;
end
ax = f(mem);
fc=1;

[ax, indx]=sort(ax);
mem=mem(indx,:);


figure
imprime(1,vxp,vyp,vzp,mem(:,1),mem(:,2),ax)


newmem(:,1)=mem(:,1);
newmem(:,2)=mem(:,2);
newmem(:,3)=ax;

diff=1;
it=0;
%% start
while abs(diff)>Converge
    it=it+1;
    lastmem(:,1)=newmem(:,1);
    lastmem(:,2)=newmem(:,2);
    lastmem(:,3)=newmem(:,3);
    
    elite1x=newmem(1:ceil(surv*N),1);
    elite1y=newmem(1:ceil(surv*N),2);
    elite1a=newmem(1:ceil(surv*N),3);
    
mem=newmem(:,1:2);

%% crossover & mutation & affinity calc
corandx=rand(round((1-crovrate)*N),1);
corandy=rand(round((1-crovrate)*N),1);
for i=1:N/2
    xbaby(i)=mem(i,1)-corandx(i)*(mem(i,1)-mem(i+1,1));%*search;
    ybaby(i)=mem(i,2)-corandy(i)*(mem(i,2)-mem(i+1,2));%*search;
    if abs(mem(i,1) - mem(i+1,1)) < 1e-4 && abs(mem(i,2) - mem(i+1,2)) < 1e-4
        if mem(i,1) < 0
            xbaby(i)=mem(i,1)+corandx(i)*search;
        elseif mem(i,1) > 0
            xbaby(i)=mem(i,1)-corandx(i)*search;
        end
        if mem(i,2) < 0
            ybaby(i)=mem(i,2)+corandy(i)*search;
        elseif mem(i,2) > 0
            ybaby(i)=mem(i,2)-corandy(i)*search;
        end
    end
end
for i=N/2+1:N
    xbaby(i)=mem(i,1)+corandx(N)*(mem(i-N/2,1)-mem(i,1));%*search;
    ybaby(i)=mem(i,2)+corandy(N)*(mem(i-N/2,2)-mem(i,2));%*search;
end


% mutation
mutpos=randi([1,N],ceil(mutrate*N),1);
xmut=varMin+(varMax-varMin).*rand(ceil(mutrate*N),1);
ymut=varMin+(varMax-varMin).*rand(ceil(mutrate*N),1);
mem(mutpos,1)=xmut;
mem(mutpos,2)=ymut;

%expectation cal for parents
axpar = f(mem);
fc=fc+1;
for i=1:N
    same=(axpar(:)==axpar(i));
    samenum=sum(same(:) == 1);
    c(i)=1-samenum/N;
    e(i)=axpar(i)/c(i);
end

% expectation cal for baby
axbaby= f([xbaby,ybaby]);
fc=fc+1;
for i=1:size(axbaby,1)
    same=(axbaby(:)==axbaby(i));
    samenum=sum(same(:) == 1);
    cbaby(i,1)=1-samenum/(size(axbaby,1));
    ebaby(i,1)=axbaby(i)/cbaby(i,1);
end

hx=[mem(:,1); xbaby];
hy=[mem(:,2); ybaby];
hax=[axpar; axbaby];
he=[e; ebaby];

[he, indx]=sort(he);
hx=hx(indx);
hy=hy(indx);
hax=hax(indx);

% affinity calc
haff=ones(N+N/2,N+N/2-1)*NaN;
cont=1;
while cont==1
for i=ceil(surv*N)+1:N+N/2-1
    for j=i+1:N+N/2
        distance=sqrt((hx(i)-hx(j))^2+(hy(i)-hy(j))^2);
        haff(j,i)=1/(1+distance);
    end
end

[temp, col]=maxk(maxk(haff,1),1);
[temp, rw]=maxk(haff,1);
row=rw(col);
haff(row,col)=NaN;
if he(row)-he(col)<0
    del=row;
else
    del=col;
end
hx(del)=NaN;
hy(del)=NaN;
he(del)=NaN;

if sum(isnan(he))==N/2
    cont=0;
end
end


% save to memory cell
[he, indxm]=sort(he); % sort descending order for maximization problem
hx=hx(indxm);
hy=hy(indxm);
mem(:,1)=hx(1:N);
mem(:,2)=hy(1:N);
hax=hax(indxm);
ax=hax(1:N);

% ax= f(mem);
% fc=fc+1;

[ax, indxsol]=sort(ax);
mem(:,1)=mem(indxsol,1);
mem(:,2)=mem(indxsol,2);

for el=1:ceil(surv*N)
    if elite1a(el)<ax(el)
        mem(el,1)=elite1x(el);
        mem(el,2)=elite1y(el);
        ax(el)=elite1a(el);
    end
end


%% update memcell
newmem=lastmem;
updated=0;
for ii=1:N
    if ax(ii) < lastmem(ii,3)
        newmem(ii,1)=mem(ii,1);
        newmem(ii,2)=mem(ii,2);
        newmem(ii,3)=ax(ii);
    end
end

imprime(1,vxp,vyp,vzp,newmem(:,1),newmem(:,2),newmem(:,3))

gen(it,:,:)=newmem;
optx(it,1)=gen(it,1,1);
opty(it,1)=gen(it,1,2);
opta(it,1)=gen(it,1,3);

%Termination Condition
if it > maxit
%     if gen(it-Termcond,1,3)==opta(it)
        disp('termination')
        break
%     end
end

if it > 5
    avgpast=sum(gen(it-5,1:ceil(0.2*N),3))/ceil(0.2*N);
    avgcurr=sum(gen(it,1:ceil(0.2*N),3))/ceil(0.2*N);
    diff=abs(avgpast-avgcurr)/avgpast*100;
end

    fprintf('%2d  f(%6.2e,%6.2e): %7.2e\n',it,optx(it,1),opty(it,1),opta(it,1))

end

toc


optval=gen(end,1,3);


%% Plot
figure(1)
imprime(1,vxp,vyp,vzp,newmem(:,1),newmem(:,2),newmem(:,3))

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
        view(0,-90) %COMMENT ME FOR 3D PLOT VIEW
        title('Minimization')
        xlabel('x_1')
        ylabel('x_2')
%         zlabel('f(x_1,x_2)') %ADD ME FOR 3D PLOT VIEW
%         plot3(x,y,fx,'k*') %ADD ME FOR 3D PLOT VIEW
        scatter(x,y,'k*') %COMMENT ME FOR 3D PLOT VIEW
        colormap jet
        drawnow
        hold off
        pause(0.1)
    end

end