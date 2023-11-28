%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Mauricelle and Vargas
% Last update: Nov 28, 2023
% Motivation: experimental data collected
% from a shaking table,  identification
% procedure for quasi-LPV of linear systems
% E-mail: avargas@utfpr.edu.br
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all, close all, clc, format long, format compact,

disp(' .... computing code for robust linear system from Shaking Table (it may take some minutes) ...')


fid = fopen('listaData.txt');
tline = fgetl(fid);
count = 1;
while ischar(tline)
    nome{count} = sprintf('%s',tline);
    tline = fgetl(fid);
    count = count+1;
end

fclose(fid);


for cx=1:max(size(nome))
    
    load(nome{cx});
    
    Ut =  out.Ut;
    
    x1real = out.timepos(:,2);
    x2real = out.Velocity;
    
    
    TS = 0.005;
    VarEPS = 1e-6;
    
    Nit = max(size(Ut(1:end)));
    Ts=TS;
    t=[0:Ts:Ts*(Nit)];
    
    P{1} = 100000*diag([VarEPS VarEPS 10^3 10^3  VarEPS 10^3]);
    Theta{1} = [0; 0; 0 ;0;0;0] ;   % initialize Theta
    Q = diag([1 1 1 1 1 1]);
    Sigma = 1;
    I = eye(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % procedure of Kalman filter to identify
    % the parameters of the linear system
    % x(k+1) = A*x(k) + B*u(k) + cte
    % We write the system above in the form of
    % \theta(k+1) = \theta(k) + w(k)
    %  x(k) = m(k)*theta(k) + e(k)
    % The last two equations feed the KF
    
    U_shifted=[0;Ut];
    x1real_shifted = [0;x1real];
    x2real_shifted = [0;x2real];
    
    Phi = zeros(2, 6, 4000);
    
    for k=1:Nit
        Phi(:,:,k) = [x1real_shifted(k) x2real_shifted(k) 0 0 U_shifted(k) 0; 0 0 x1real_shifted(k) x2real_shifted(k) 0 U_shifted(k)];
    end
    
    for k=2:Nit
        y{k} = [x1real(k); x2real(k)];
        K{k} = P{k-1}*Phi(:,:,k)'*(Phi(:,:,k)*P{k-1}*Phi(:,:,k)' + Sigma*I)^(-1);
        Theta{k} = Theta{k-1} + K{k}*(y{k} - Phi(:,:,k)*Theta{k-1});
        P{k} = P{k-1} - K{k}*Phi(:,:,k)*P{k-1} + Q;
        P{k} = (P{k}+P{k}')/2;  %use this line to increase stability
    end
    
    vecx1=[];
    vecx2=[];
    xsim = [x1real(1) x2real(1)]';
    % xsim{1}= [x1real(1) x2real(1)]';
    
    for k=1:Nit
        
        vecx1=[vecx1  xsim(1)];
        vecx2=[vecx2  xsim(2)];
        
        % WRITE A CODE THAT WILL SIMULATE
        % YOUR ESTIMATED MATRICES A{k}, B{k}
        
        Theta1(:,:,k) = Theta{k}; % converter cell para double
        
        A = [Theta1(1,1,k) Theta1(2,1,k)
            Theta1(3,1,k) Theta1(4,1,k)];
        B  = [Theta1(5,1,k) Theta1(6,1,k)]';
        
        xsim =[Theta1(1,1,k)*x1real(k)+ Theta1(2,1,k)*x2real(k) +  Theta1(5,1,k)*Ut(k);  Theta1(3,1,k)*x1real(k) +  Theta1(4,1,k)*x2real(k) + Theta1(6,1,k)*Ut(k)];
        
    end
    
    a11=[];
    a12=[];
    a21=[];
    a22=[];
    b1=[];
    b2=[];
    for k=1:Nit
        a11 = [a11; Theta1(1,1,k)];
        a12 = [a12; Theta1(2,1,k)];
        a21 = [a21; Theta1(3,1,k)];
        a22 = [a22; Theta1(4,1,k)];
        b1 = [b1; Theta1(5,1,k)];
        b2 = [b2; Theta1(6,1,k)];
    end
    
    % RMSE
    sum1 = 0;
    sum2 = 0;
    for k = 1:4001
        sum1 = sum1 + (x1real(k) - vecx1(k))^2;
        sum2 = sum2 + (x2real(k) - vecx2(k))^2;
    end
    RMSE1 = sqrt(sum1/k);
    RMSE2 = sqrt(sum2/k);
    
    %% Création d'un cell array pour stocker les matrices A
    num_matrices = 4001;
    matrices_A = cell(1, num_matrices);
    matrices_B = cell(1, num_matrices);
    
    for k = 1:num_matrices
        
        aij_values = [a11(k) a12(k); a21(k) a22(k)];
        bi_values = [b1(k);b2(k)];
        
        A = reshape(aij_values, 2, 2);
        B = reshape(bi_values, 2, 1);
        
        matrices_A{k} = A;
        matrices_B{k} = B;
    end
    
    savefile = sprintf('rough_matrices_%0.3i.mat',cx)
    save(savefile, 'matrices_A', 'matrices_B','-v7');
    
    
    close all,
    
    figure(1)
    subplot(2,1,1)
    hold on
    plot(Ts*[1:max(size(x1real))],x1real,'r-.','LineWidth',2);
    plot(Ts*[1:max(size(x1real))],vecx1,'k','LineWidth',1)
    hold off
    grid
    legend('real','estimated');
    xlabel('seconds'),ylabel('Displacement(x1)')
    
    subplot(2,1,2)
    hold on
    plot(Ts*[1:max(size(x2real))],x2real,'r-.','LineWidth',2);
    plot(Ts*[1:max(size(x2real))],vecx2,'k','LineWidth',1)
    hold off
    grid
    legend('real','estimated');
    xlabel('seconds'),ylabel('Velocity(x2)')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    
    subplot(4,2,1)
    plot(Ts*[1:max(size( Ut ))], Ut ,'b','LineWidth',2);
    ylabel('Ut'),
    grid
    
    subplot(4,2,2)
    plot(Ts*[1:max(size(a11))],a11,'k','LineWidth',2);
    grid, ylabel('a11'),
    
    subplot(4,2,3)
    plot(Ts*[1:max(size( a12 ))],a12,'k','LineWidth',2);
    grid, ylabel('a12'),
    
    subplot(4,2,4)
    plot(Ts*[1:max(size( a21 ))],a21,'k','LineWidth',2);
    grid, ylabel('a21'),
    
    subplot(4,2,5)
    plot(Ts*[1:max(size( a22 ))],a22,'k','LineWidth',2);
    grid, ylabel('a22'),
    
    subplot(4,2,6)
    plot(Ts*[1:max(size( b1 ))],b1,'k','LineWidth',2);
    grid, ylabel('b1'),
    
    subplot(4,2,7)
    plot(Ts*[1:max(size( b2 ))],b2,'k','LineWidth',2);
    grid, ylabel('b2');
    
end

