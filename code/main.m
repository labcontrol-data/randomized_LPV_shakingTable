%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Mauricelle and Vargas
% Last update: Jan 31, 2024
% Motivation: experimental data collected
% from a shaking table. Procedure that computes 
% stability according to `randomized approach'
% E-mail: avargas@utfpr.edu.br
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all, close all, clc, format long, format compact,

disp(' .... procedure to compute randomized stability test (it may take some hours) ...')

reset(RandStream.getGlobalStream,sum(100*clock));  % ensures random seed to Matlab

Nit=3;  %horizon or number of steps

Rep = 5;

text_file = sprintf('vertices_final.mat');
load(text_file);

A = A_vertices;
B = B_vertices;

N=max(size(A)); % number of vertices
[n,l]=size(B{1}); %dimension of matrix B

vecFeas = [];

varEps = 1e-4;
c1 = 1;
c2 = 1;
xi = ((1 - varEps)^N)*c1/c2  % xi has to be less than one

for m=1:Rep
    m
   
    warning('off','YALMIP:strict')
    %##########################################################################
    %#  LMIS modificadas ....
    %##########################################################################
    LMI=set([]);
    for i=1:N
        P{i}=sdpvar(n,n,'sy');
        LMI = [LMI, P{i}>= 0]; ;
    end
   
    K =[-0.2   2];
   
    for k=1:Nit
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
   
    vecFeas = [];
    for k=1:Nit
        if (k==Nit)
            LMI = [LMI, P_alpha{1} < c1*eye(n,n)];
            LMI = [LMI, P_alpha{Nit} > c2*eye(n,n)];
        else
            A_cl = A_alpha{k} + B_alpha{k}*K;
            LMI = [LMI, A_cl'*P_alpha{k+1}*A_cl - (1 - varEps)*P_alpha{k}  < 0 ];
        end
    end
   
    ops=sdpsettings('solver','mosek','verbose',1);
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
