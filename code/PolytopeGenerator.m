function [Av,Bv] = PolytopeGenerator(A,B)

disp('.... procedure may take many hours or days to complete  ....')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the 'M' matrices. Both A(alpha) and B(alpha)
% have the same index 'alpha'. This is the reason
% A(alpha) and B(alpha) must be together to create 'M'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=[];
for k=1:max(size(A))
    M{k} = [A{k}  B{k}];
end
    
r = 9;  % it can be any number
Delta = {M{r+1},M{r}};

for j=1:max(size(A))
    j
    n = length(Delta);
    LMI = [];
    for i=1:n
        alpha{i} = sdpvar(1,1,'full');
        LMI = [LMI, alpha{i} >= 0];
    end

    M_alpha = 0; sx = 0;
    for i=1:n
        M_alpha = M_alpha + alpha{i}*Delta{i};
        sx = sx + alpha{i};
    end
    LMI = [LMI, sx == 1 ];
    LMI = [LMI, M{j} == M_alpha ];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % objetivo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    ops=sdpsettings('solver','mosek','verbose',0);
    quiz = solvesdp(LMI,[],ops);
    p = min(checkset(LMI));
    disp(' ');
    
    if (quiz.problem == 0)   % other condition (more imprecise)
                             % is to consider (p>=0)
        disp('YES:Feasible')
        text = [];
        for i=1:n
            text = [text sprintf('alpha(%d)=%0.5g; ',i,double(alpha{i}))];
        end
        disp(text)
    else
        disp('NO:Unfeasible')
        Delta = [Delta, {M{j}}];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the set of vertices 'Av' and 'Bv' 
% for the polytopes A(alpha) and B(alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[la,ca] = size(A{1});
[lb,cb] = size(B{1});
Av=[]; Bv=[];
for k=1:max(size(Delta))
    Av{k} = Delta{k}(1:la,1:ca);
    Bv{k} = Delta{k}(1:la,ca+1:ca+cb);
end

