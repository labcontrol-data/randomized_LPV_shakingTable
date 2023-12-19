%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Mauricelle and Vargas
% Last update: Dec 18, 2023
% Motivation: experimental data collected
% from a shaking table,  identification
% procedure for quasi-LPV of linear systems
% E-mail: avargas@utfpr.edu.br
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all, close all, clc, format long, format compact,

disp(' .... identification for Shaking Table (it may take some minutes) ...')


fid = fopen('listaData.txt');
tline = fgetl(fid);
count = 1;
while ischar(tline)
    nome{count} = sprintf('%s',tline);
    tline = fgetl(fid);
    count = count+1;
end

fclose(fid);

factor = [10 1];     %values obtained by trial and error 

n = 1;
cx = 1;
while cx<=max(size(nome))
    
    load(nome{cx});
    
    Ut =  out.Ut;
    u = Ut;
    
    x1real = out.timepos(:,2)/100;
    x2real = out.Velocity;
    
    TS = 0.005;
    
    Ts=TS;
    
    VarEPS = 1;
    
    Nit = max(size(u));
    Ts=TS;
    t=[0:Ts:Ts*Nit];
    
    P{1} = 1*eye(7,7);
    Delta{1} = [0 ; 0 ; 0; 0; 0; 0; 0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % procedure of Kalman filter to identify
    % the parameters of the linear system
    % x(k+1) = A*x(k) + B*u(k) + cte
    % We write the system above in the form of
    % \theta(k+1) = \theta(k) + w(k)
    %  x(k) = m(k)*theta(k) + e(k)
    % The last two equations feed the KF
    for k=2:Nit
        m{k-1} = [x1real(k-1) x2real(k-1)  0   0  u(k-1)  0 1;
            0  0  x1real(k-1) x2real(k-1) 0 u(k-1) 0];
        K{k} = P{k-1}*m{k-1}'*inv( factor(2)*eye(2,2) + m{k-1}*P{k-1}*m{k-1}' ) ;
        Delta{k} = Delta{k-1} + K{k}*( [x1real(k) x2real(k)]' - m{k-1}*Delta{k-1} );
        P{k} = factor(1)*diag([1, 1 ,1,1,1,1,1])  + P{k-1} - K{k}*m{k-1}*P{k-1};
        P{k} = (P{k}+P{k}')/2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    heap_a11=[]; heap_a12=[];
    heap_a21=[]; heap_a22=[];
    heap_b1=[]; heap_b2=[]; heap_h=[];
    vecx1=[];vecx2=[];
    xsim{1}= [x1real(1) x2real(1)]';
    for k=1:Nit
        vecx1 = [vecx1 xsim{k}(1)];
        vecx2 = [vecx2 xsim{k}(2)];
        A = [Delta{k}(1:2,:)';
            Delta{k}(3:4,:)'];
        B = [Delta{k}(5:6,:)];
        H = [Delta{k}(7,:); 0];
        xsim{k+1} = A*xsim{k} + B*u(k) + H*1;
        
        heap_a11 = [heap_a11  A(1,1)];
        heap_a12 = [heap_a12  A(1,2)];
        heap_a21 = [heap_a21  A(2,1)];
        heap_a22 = [heap_a22  A(2,2)];
        heap_b1 = [heap_b1 B(1,1)];
        heap_b2 = [heap_b2 B(2,1)];
        heap_h = [heap_h H(1,1)];
    end
    
    
    for k=1:Nit
        A_d{k} = [heap_a11(k) heap_a12(k);
            heap_a21(k) heap_a22(k)];
        B_d{k} = [heap_b1(k); heap_b2(k)];
        H_d{k} = [heap_h(k); 0];
    end    
    
    error = norm(vecx1-x1real);
    cx
    if (error>1000)
        disp(nome{cx});
        factor = [0.1*rand 10];   %values obtained by trial and error
    else    
        factor = [10 1];       %values obtained by trial and error
               
        figure(n)
        subplot(2,1,1)
        hold on
        plot(Ts*[1:max(size(x1real))],x1real,'r-.','LineWidth',2);
        plot(Ts*[1:max(size(x1real))],vecx1,'k','LineWidth',1)
        hold off
        grid
        legend('real','estimated');
        xlabel('seconds'),ylabel('x1')
        
        subplot(2,1,2)
        hold on
        plot(Ts*[1:max(size(x2real))],x2real,'r-.','LineWidth',2);
        plot(Ts*[1:max(size(x2real))],vecx2,'k','LineWidth',1)
        hold off
        grid
        legend('real','estimated');
        xlabel('seconds'),ylabel('x2')
        
        cx = cx+1;
        
        matrices_A = A_d;
        matrices_B = B_d;
        matrices_H = H_d;
    
        savefile = sprintf('rough_matrices_%0.3i.mat',n)
        save(savefile, 'matrices_A', 'matrices_B', 'matrices_H','-v7');
        
        n = n+1;
        
    end
    
end
