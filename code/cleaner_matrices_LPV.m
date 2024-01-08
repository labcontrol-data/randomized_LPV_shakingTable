%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Mauricelle and Vargas
% Last update: Jan 8, 2024
% Motivation: experimental data collected
% from a shaking table. Procedure that "clears out"
% matrices that are too close with each other.
% E-mail: avargas@utfpr.edu.br
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all, clc, format short, format compact,


disp(' .... procedure for clean matrices: Shaking Table (it may take some minutes) ...')


varEps = 0.01; % tolerance, usually small, like 0.01 or 0.001

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
    disp(cx)        
    
    text_file = sprintf('rough_matrices_%0.3i.mat',cx);
    load(text_file);
    
    A_rough = matrices_A;
    B_rough = matrices_B;
    
    init = 50; % this number can be any
    A_po = []; B_po = []; 
    A_po{1} = A_rough{init};
    B_po{1} = B_rough{init};
    
    for j=1:max(size(A_rough))
        vecFlag = [];   % flag indicates that two numbers are close (0)
                        % or away (1) one from another
        count = 1;
        while (count<=max(size(A_po)))
            if norm([A_rough{j} B_rough{j}] -   [A_po{count} B_po{count}] ,'fro')<varEps
                vecFlag = [vecFlag 0];
            else
                vecFlag = [vecFlag 1];
            end
            count = count+1;
        end
        if (sum(vecFlag)==max(size(vecFlag)))
            A_po = [A_po, A_rough{j}];
            B_po = [B_po, B_rough{j}];
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % procedure to 'clean' the rough terms from `(A_rough,B_rough)'
    % Note that `(A_clean,B_clean)' and `(A_rough,B_rough)' have the 
    % same number of elements. Every element of `(A_clean,B_clean)' is 
    % at least `varEps' far away from each other
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_clean{1}  = A_rough{1};
    B_clean{1}  = B_rough{1};
    for j=2:max(size(A_rough))
        A_clean{j} = A_rough{j};
        B_clean{j} = B_rough{j};
        count=1;
        while (count<=max(size(A_po)))
            if norm([A_rough{j} B_rough{j}] - [A_po{count} B_po{count}],'fro')<varEps
                A_clean{j} = A_po{count};
                B_clean{j} = B_po{count};
                count = max(size(A_po));
            end
            count = count+1;
        end
    end
    
    matrices_A = A_clean;
    matrices_B = B_clean;
        
    savefile = sprintf('clean_matrices_%0.3i.mat',cx);
    save(savefile, 'matrices_A', 'matrices_B','A_po','B_po','-v7');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % procedure to check how the 'clean procedure' affects
    % the system representation --- procedure computes error
    % to see if it justifies a Gaussian-like curve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load(nome{cx});
    u =  out.Ut;
    x1real = out.timepos(:,2)/100;
    x2real = out.Velocity;
    TS = 0.005;
    Ts=TS;
       
    text_file = sprintf('clean_matrices_%0.3i.mat',cx);
    load(text_file);
    A_clean = matrices_A;
    B_clean = matrices_B;
        
    text_file = sprintf('rough_matrices_%0.3i.mat',cx);
    load(text_file);
    A_rough = matrices_A;
    B_rough = matrices_B;
    H_rough = matrices_H;
        
    vecx1=[];vecx2=[]; error=[];
    xsim{1}= [x1real(1) x2real(1)]';
    for k=1:max(size(x1real))
        vecx1 = [vecx1 xsim{k}(1)];
        vecx2 = [vecx2 xsim{k}(2)];
        A = A_rough{k};
        B = B_rough{k};
        H = H_rough{k};
        xsim{k+1} = A*xsim{k} + B*u(k) + H*1;
        error{k} = (A_rough{k} - A_clean{k})*xsim{k} +...
                   (B_rough{k} - B_clean{k})*u(k) + ...
                   (H_rough{k}); % - H_clean{k});
    end
    heap_error{cx} = error;
     
%     figure(cx)
%     subplot(2,1,1)
%     hold on
%     plot(Ts*[1:max(size(x1real))],x1real,'r-.','LineWidth',2);
%     plot(Ts*[1:max(size(x1real))],vecx1,'k','LineWidth',1)
%     hold off
%     grid
%     legend('real','estimated');
%     xlabel('seconds'),ylabel('x1')
%     
%     subplot(2,1,2)
%     hold on
%     plot(Ts*[1:max(size(x2real))],x2real,'r-.','LineWidth',2);
%     plot(Ts*[1:max(size(x2real))],vecx2,'k','LineWidth',1)
%     hold off
%     grid
%     legend('real','estimated');
%     xlabel('seconds'),ylabel('x2')
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% procedure that shows the error like a Gaussian curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vecE1=[];

vecE2=[];
for cx=1:max(size(nome))
    error = heap_error{cx};
    for j=1:max(size(error))
        vecE1 = [vecE1 error{j}(1)];
        vecE2 = [vecE2 error{j}(2)];
    end
end


figure(10*cx)
subplot(2,1,1)
plot(vecE1,'b','LineWidth',1);
grid, legend('error');
xlabel('sample'),ylabel('error')

subplot(2,1,2)
plot(vecE2,'b','LineWidth',1);
grid, legend('error');
xlabel('sample'),ylabel('error')


figure(10*cx+1)
histogram(vecE1)

figure(10*cx+2)
histogram(vecE2)
