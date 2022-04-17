                %  -------------------------------  %
                %  Computational E/M                %
                %   Angelitsi Sotiria, AEM:4366     %
                %   Set 2                           %
                %  -------------------------------  %
                
clear; clc; clf;
c=1;    % speed of light
phi=0:0.25:90; % propagation angle
N=361; % diamerisi tou axona twn goniwn
Nlamda=[5.0 10.0 20.0]; % plithos shmeiwn ana m.k.
lamda0=1.0; % m.k.
S=0.5; % paragontas gia amesi sygrish ths FDTD me to thewrhtiko apotelesma

%% Euresi analogias vp/c(phi) gia omoiomorfo plegma (dx=dy=d)

for i=1:3 %epanalhptiki methodos gia tis 3 times tou Nlamda
    d(i)=lamda0./Nlamda(i); %Dx=Dy=D
    dx(i)=d(i);
    dy(i)=d(i);
    for j=1:N  % epanalhptikothta gia kathe gwnia pou tou ypodiksame
        A=dx(i)*cos(phi(j)*pi/180)/2; % oi oroi gia th sxesi tou k apo Taflove
        B=dy(i)*sin(phi(j)*pi/180)/2;
        C=sin(pi*S/Nlamda(i))*sin(pi*S/Nlamda(i))/(S*S);
        ki=2*pi; % arxiki timi tou k
        for m=1:3 % Epanalhptikh methodos Newton-Raphson. 3 epanalipseis einai arketes.
          kiplus1=ki-(sin(A*ki)^2+sin(B*ki)^2-C)/(A*sin(2*A*ki)+B*sin(2*B*ki));
        end
        ki=kiplus1; % swap gia euresi ths epomenhs timhs
        vpratio(j)=2*pi/kiplus1; % Ypologismos vp/c synarthsei tou phi
    end
    f1=figure(1)
    f1.Name=('Uniform Grid');
    set(f1,'NumberTitle', 'off');
    plot(phi,vpratio(1,:),'-.')
    title('Numerical Dispersion in 2D implementation of FDTD by Newton΄s iterative method')
    xlabel('Wave angle, φ (degrees)')
    ylabel('Normalized phase velocity, vp/c')
    legend({'lo/5','lo/10','lo/20'},'location','south')
    hold on
end

%% Euresi analogias vp/c(phi) gia anomoiomorfo plegma (dx=2dy=d)
for i=1:2 %epanalhptiki methodos gia tis 3 times tou Nlamda
    d2=1.0/Nlamda(i); %Dx=2Dy=D
    dx2=d2;
    dy2=d2/2.;
    for j=1:N  % epanalhptikothta gia kathe gwnia pou tou ypodiksame
        ki2=2*pi;  % arxiki timi tou k
        for m=1:3 % Epanalhptikh methodos Newton-Raphson. 3 epanalipseis einai arketes.
            % Ypologismos tou k apo th sxesh stis diafaneies tou Gedney
            kiplus12=ki2-((sin(ki2*cos(phi(j)*pi/180)*dx2/2))^2/dx2^2+ ...
                (sin(ki2*sin(phi(j)*pi/180)*dy2/2))^2/dy2^2 ...
                -(1/(S^2*dy2^2))*sin(pi*S*dy2)^2)/ ...
                ((cos(phi(j)*pi/180)/(2*dx2))*sin(ki2*cos(phi(j)*pi/180)*dx2)+ ...
                sin(phi(j)*pi/180)*sin(ki2*sin(phi(j)*pi/180)*dy2)/(2*dy2));
        end
        ki2=kiplus12; % swap gia euresi ths epomenhs timhs
        vpratio2(j)=2*pi./kiplus12; % Ypologismos vp/c synarthsei tou phi
    end
    f2=figure(2)
    f2.Name=('Non-Uniform Grid');
    set(f2,'NumberTitle', 'off');
    plot(phi,vpratio2(1,:),'-.')
    title('Numerical Dispersion in 2D implementation of FDTD by Newton΄s iterative method')
    xlabel('Wave angle, φ (degrees)')
    ylabel('Normalized phase velocity, vp/c')
    legend({'lo/5','lo/10'},'Location','east')
    hold on
end
