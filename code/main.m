%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Mauricelle and Vargas
% Last update: May 23, 2024
% Motivation: experimental data collected
% from a shaking table. Procedure that xxxxxxxxx
% xxxxxxxxxxxxxxxxxxxx.
% E-mail: avargas@utfpr.edu.br
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all, close all, clc, format long, format compact,

disp(' .... procedure to compute randomized stability test (it may take some hours) ...')

reset(RandStream.getGlobalStream,sum(100*clock));  % ensures random seed to Matlab

T=8;  %horizon or number of steps

Rep = 10;  %number of repetitions for the randomized approach

text_file = sprintf('vertices_final.mat');
load(text_file);

A = A_vertices;
B = B_vertices;
C = [1  0];

N=max(size(A));   % number of vertices
[n,l]=size(B{1}); %dimension of matrix B

vecFeas = [];

c0 = 1e-4;
c1 = 0.998;
c2 = 0.999;
xi = ((1 - c0)^T)*c1/c2  % xi has to be less than one

K1 = 0.5;
K2 = 0.1; 
K3 = 60; 

for m=1:Rep
    m
   
    warning('off','YALMIP:strict')
    %##########################################################################
    %#  main LMIS 
    %##########################################################################
    LMI=set([]);
    for i=1:N
        P{i}=sdpvar(n+1,n+1,'sy');
        LMI = [LMI, P{i}>= 0]; ;
    end
   
    for k=1:T
        A_alpha{k} = 0;
        B_alpha{k} = 0;
        P_alpha{k} = 0;
       
        aleat = rand(1,N);
        alpha = aleat/sum(aleat);
       
        for i=1:N
            A_alpha{k} = A_alpha{k} + alpha(i)*A{i};
            B_alpha{k} = B_alpha{k} + alpha(i)*B{i};
            P_alpha{k} = P_alpha{k} + alpha(i)*P{i};
        end
    end
   
    for k=1:T
        if (k==T)
            LMI = [LMI, P_alpha{1} < c1*eye(n+1,n+1)];
            LMI = [LMI, P_alpha{T} > c2*eye(n+1,n+1)];
        else
            % A_cl = (A_alpha{k} + B_alpha{k}*K);
            A_cl = [A_alpha{k} + B_alpha{k}*[K1 K2], B_alpha{k}*K3; -C, eye(l,l)];
            LMI = [LMI, A_cl'*P_alpha{k+1}*A_cl - (1 - c0)*P_alpha{k}  < 0 ];
        end
    end
   
    ops=sdpsettings('solver','mosek','verbose',0);
    quiz = solvesdp(LMI,[],ops);
    p=min(checkset(LMI));
   
    %capturing the solutions (if they exist)
   
    if ((p > -1e-12)&&(quiz.problem==0))
        disp('LMI feasible')
        vecFeas = [vecFeas 0];
    else
        disp('LMI unfeasible')
        vecFeas = [vecFeas 1];
    end
   
end

if (sum(vecFeas)==0)
    disp('---> All LMIs are feasible (good for us) <---')
else
    disp('One LMI was infeasible. Sorry.')
end
