function [ineq_iis,eq_iis,ncalls] = deletionfilter(Aineq,bineq,Aeq,beq,lb,ub)
ncalls = 0;
[mi,n] = size(Aineq); % Number of variables and linear inequality constraints
f = zeros(1,n);
me = size(Aeq,1); % Number of linear equality constraints
opts = optimoptions("linprog","Algorithm","dual-simplex","Display","none");

ineq_iis = true(mi,1); % Start with all inequalities in the problem
eq_iis = true(me,1); % Start with all equalities in the problem

for i=1:mi
    ineq_iis(i) = 0; % Remove inequality i
    [~,~,exitflag] = linprog(f,Aineq(ineq_iis,:),bineq(ineq_iis),...
                                Aeq,beq,lb,ub,opts);
    ncalls = ncalls + 1;
    if exitflag == 1 % If now feasible
        ineq_iis(i) = 1; % Return i to the problem
    end
end
for i=1:me
    eq_iis(i) = 0; % Remove equality i
    [~,~,exitflag] = linprog(f,Aineq,bineq,...
                                Aeq(eq_iis,:),beq(eq_iis),lb,ub,opts);
    ncalls = ncalls + 1;
    if exitflag == 1 % If now feasible
        eq_iis(i) = 1; % Return i to the problem
    end
end
end