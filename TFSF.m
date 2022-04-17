clear; clc; clf;

N=200;                  % Number of Cells
c=3e+8;                 % Speed of light

% Free space
mi0=4.0*pi*1.0e-7;      % Permeability
epsilon0=1.0/(c*c*mi0); % Permittivity

% Material parameters
epsilon=1.0*ones(1,N);         % Relative electric permittivity
sigmaE=0.0*ones(1,N);          % Electric conductivity
epsilon(100:160)=9;            % Relative electric permittivity (PEC)
updatesigmaE=[0.01 0.5 1]; % Electric conductivity (PEC)                 
mi=1.0*ones(1,N);
sigmaH=0.0*ones(1,N);

%Source
updatefreq=[1.0e+9 3.0e+9 5.0e+9];            % Frequency

%% Update Coefficients in Yee 1-D scheme
for count1=1:6 % From 1-3 sigmaE changes. From 4-6 frequency changes
    if (count1<=3)
        sigmaE(100:160)=updatesigmaE(count1); 
        freq=1.0e+9;
    else
        sigmaE(100:160)=0.01;
        count2=count1-3;
        freq=updatefreq(count2);
    end
    
    %Source
    S=1.0;                  % cdt/dx (magic step)
    lambda=c/freq;          % Wavelength
    dx=lambda/20.0;         % 20 cells per wavelength
    dt=S*dx/c;              % Time step
    
    Ex0=1.0;                % Initialization
    nmax=3*round(N*S);      % Total time for computing and animation
    source(1:nmax)=Ex0*sin(2*pi*freq*(1:nmax)*dt); % Source
    
    % Initial Conditions
    Ex(1:N)=0.0;
    Hy(1:N-1)=0;
    emaxE(1:N)=0.0;
    emaxH(1:N-1)=0.0;
    EmaxTot(N,nmax)=0.0;

    for i=1:N
        eprop=sigmaE(i)*dt/(2.0*epsilon0*epsilon(i));
        ca(i)=(1.0-eprop)/(1.0+eprop);
        cb(i)=dt/(epsilon0*epsilon(i)*dx)/(1.0+eprop);
        hprop=sigmaH(i)*dt/(2.0*mi0*mi(i));
        da(i)=(1.0-hprop)/(1.0+hprop);
        db(i)=dt/(mi0*mi(i)*dx)/(1.0+hprop);
    end    

%% Time Stepping - Yee Algorithm - 1D FDTD Method

    for j=1:nmax
        Ex(1)=source(j);
    
        Hy(1:N-1)=da(1:N-1).*Hy(1:N-1)+db(1:N-1).*(Ex(2:N)-Ex(1:N-1));

        Ex(2:N-1)=ca(2:N-1).*Ex(2:N-1)+cb(2:N-1).*(Hy(2:N-1)-Hy(1:N-2));
        E(N)=0;
    
        emaxE=max(Ex,emaxE);   % Megista plath gia kathe xroniki stigmi
        emaxH=max(Hy,emaxH);
        Ex(1)=Ex(2);
        Ex(N)=Ex(N-1);
        
        EmaxTot(:,j)=emaxE;    % Ola ta megista plath gia tous 3 gyrous
    
        if (count1==1) % Gia na trexei mono to 1o stigmiotypo.AN AFAIRAITHEI FAINONTAI KAI TA 6!
            f1=figure(1);   % Animation Ex,Hy (Grid Coordinates)
            f1.Name=['Angelitsi Sotiria - Animation:',num2str(count1)];
            set(f1,'NumberTitle', 'off');
            subplot(2,1,1),plot(1:N,Ex,'r',1:N,emaxE(1:N),'g'),axis([1 N -3 3]);
            ylabel('E_x');
            title(['Time=',num2str(round(j*dt*freq,1)),'nsec, ','σ_E=',...
                num2str(sigmaE(100)),', f=',num2str(freq/10^9),'GHz'],'fontsize',13);
    
            subplot(2,1,2),plot(1:N-1,Hy,'b',1:N-1,emaxH(1:N-1),'g'),axis([1 N-1 -8e-3 8e-3]);
            ylabel('H_y'); xlabel('Grid Coordinates')
            pause(0.05)
        end
    end

%% Computing of SkinDepth (δ)
    flag=0;
    for j=1:nmax
        for i=1:N
            if (flag==0 && EmaxTot(i,j)<1/exp(1))
                timeSD=j*dt*freq;
                delta(count1)=sqrt(2/(2*pi*freq*timeSD*mi0*sigmaE(100)));
                flag=1;
            end
        end
    end
    
end    

%% Plotting

f2=figure(2);
f2.Name=['Angelitsi Sotiria - Plot:',num2str(1)];
set(f2,'NumberTitle', 'off');
plot(updatesigmaE,delta(1:3),'*')
hold on
fplot(@(s) sqrt(1/(pi*1*1.0e+9*mi0*s)),[0.01 1])
title(['Skin Depth - ','f=',num2str(1),'GHz'],'fontsize',13)
legend('delta')
xlabel('σ_E')
ylabel('δ (μm)')

f3=figure(3);
f3.Name=['Angelitsi Sotiria - Plot:',num2str(2)];
set(f3,'NumberTitle', 'off');
plot(updatefreq,delta(4:6),'*')
hold on
fplot(@(f) sqrt(2/(2*pi*f*mi0*0.01)),[1*10^9 5*10^9])
title(['Skin Depth - ','σ_Ε=',num2str(0.01)],'fontsize',13)
legend('delta')
xlabel('Frequancy (Hz)')
ylabel('δ (μm)')
