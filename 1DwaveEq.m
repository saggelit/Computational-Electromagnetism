clear; clc; clf;
c=1;    %speed of light
Lx=10.0;    % space-steps
N=100;      % 100 shmeia gia diamerismo tou Lx
dx=Lx/(N-1);    % diamerismos twn space-steps se 100 isapexonta shmeia
time=60;    % xroniko orio pou maw endiaferei
s=[0.9 1.0 1.1];    % pinakas timwn gia ton kathorismo tou xronikou bhmatos

for i=1:3
    
    dt=s(i)*dx/c; % xroniko vhma
    
    L=zeros(N,N);           % arxikopoihsh tou tridiagwniou pinaka L
	L(  1:1+N:N*N  )=2;     % kathorismos twn timwn sth diagwnio
	L(N+1:1+N:N*N  )=-1;    % sthn apo panw grammh kai
	L(  2:1+N:N*N-N)=-1;    % sthn apo katw
    I=eye(N);               % tautotikos pinakas I me monadiaia diagwnio
    
    %% Explicit Method
    
    u_ex_nminus1=zeros(N,1);    % arxikopoihsh tou un-1
    u_ex_nminus1(2:11,1)=1;     % arxikes synthhkes
    u_ex_n=zeros(N,1);          % arxikopoihsh tou un-0
    u_ex_n(3:12,1)=1;           % arxikes synthhkes
    u_ex_nplus1=zeros(N,1);     % arxikopoihsh tou un+1
    
    A_ex=I-(c*dt)^2*L/(2*dx^2); % kathorismos tou pinaka A
    
    t_step_ex=1;
    for t_ex=1:time
        for k=2:N-1
            u_ex_nplus1=2*A_ex*u_ex_n-I*u_ex_nminus1; % Euresh lushs
        end
        if t_ex==t_step_ex*10+10 %xroniko vhma gia na krataw ta stigmiotupa 20,30...
            u_time_ex(:,t_step_ex)=u_ex_nplus1;
            t_step_ex=t_step_ex+1;
        end
        
        u_ex_nminus1=u_ex_n;   % kanw swap gia na proxwrisw xronika
        u_ex_n=u_ex_nplus1;    % n-1 <- n0   &&& n0 <- n+1
    end
    
    %% Implicit Method
    u_im_nminus1=zeros(N,1);    % arxikopoihsh tou un-1
    u_im_nminus1(2:11,1)=1;     % arxikes synthhkes
    u_im_n=zeros(N,1);          % arxikopoihsh tou un-0
    u_im_n(3:12,1)=1;           % arxikes synthhkes
    u_im_nplus1=zeros(N,1);     % arxikopoihsh tou un+1
    
    b=0.25;
    A_im=(b*L)+(dx^2/((c*dt)^2))*I;     % Kathorismos twn pinakwn A kai B
    B_im=((2*b-1)/2)*L+(dx^2/(c*dt)^2)*I;
    % C_im=A_im; %  O pinakas C de xreiasthke kathws einai idios me ton A
    
    t_step_im=1;
    for t_im=1:time
        for l=2:N-1
            u_im_nplus1=2*(A_im\B_im)*u_im_n-u_im_nminus1; % Euresh lushs
        end
        if t_im==t_step_im*10+10 %xroniko vhma gia na krataw ta stigmiotupa 20,30...
            u_time_im(:,t_step_im)=u_im_nplus1;
            t_step_im=t_step_im+1;
        end
        
        u_im_nminus1=u_im_n;      % kanw swap gia na proxwrisw xronika
        u_im_n=u_im_nplus1;       % n-1 <- n0   &&& n0 <- n+1
    end
    
    
     
   %%  Plots 
   % Dimiourgia grafimatwn se epanalhptiki diadikasia gia kathe xroniko
   % vima dt
    x=linspace(0,Lx,N);
    f1=figure(i)
    str1=s(i);
    f1.Name=['Explicit Method for Δt=',num2str(str1),'Δx/c'];
    set(f1,'color',[0.9022 0.9604 0.9703],'NumberTitle', 'off');
    for j=1:5
        subplot(5,1,j)
        plot(x,u_time_ex(:,j),'LineWidth',1.5)
        %title(['Explicit Method for Δt=',num2str(str1),'Δx/c']);
        tt1=(j+1)*10;
        legend(['t=',num2str(tt1),'sec'])
    end
    
    f2=figure(i+3)
    str2=s(i);
    f2.Name=(['Implicit Method for Δt=',num2str(str2),'Δx/c']);
    set(f2,'color',[0.9422 0.9004 0.9703],'NumberTitle', 'off');
    for j=1:5
        subplot(5,1,j)
        plot(x,u_time_im(:,j))
        str2=s(i);
        %title(['Implicit Method for Δt=',num2str(str2),'Δx/c']);
        tt2=(j+1)*10;
        legend(['t=',num2str(tt2),'sec'])
    end
end
