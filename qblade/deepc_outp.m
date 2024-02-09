function uOpt = deepc_outp(data,rf,controlParams,method,ivFlag,previewFlag)
% deepc(data,rf,controlParams,method,ivFlag,previewFlag) Function used to
% compute the DeePC optimal input.
%
% Input arguments:
%-----------------
% data          : structure containing the Hankel matrices
%                 Up,Yp,Uf,Yf,Wp,Wf (past inputs, past outputs, future
%                 inputs, future outputs, past disturbance, and future
%                 disturbance respectively), uini,yini,wini,wf (past
%                 input, output and distubance trajectories and disturbance
%                 preview respectively). Wp,Wf,wini,wf are not used if
%                 previewFlag is set to 0.
% rf            : vector of reference values for the next f steps.
% controlParams : structure containing the DeePC output weight Q and
%                 input weight R, input lower bound lbu, and input upper
%                 bound ubu, the length of the data set N, past window size
%                 p, and future window size f.
% method        : Optimization method (1 - quadprog, 2 - casadi+nlp).
% ivFlag        : Flag that indicates whether or not instrumental variables
%                 are used (1 - IV used, 0 - IV not used)
% previewFlag   : Flag that indicates whether DeePC uses preview
%                 information of the incoming disturbance or not
%                 (1 - preview used, 0 - preview not used).
%
% Output arguments:
%------------------
% uOpt          : Optimal control input.
%==========================================================================
Up = data.Up;
Yp = data.Yp;
Uf = data.Uf;
Yf = data.Yf;
uini = data.uini;
yini = data.yini;

% Input/output samples within past/future windows
pIn = size(Up,1);
pOut = size(Yp,1);
fIn = size(Uf,1);
fOut = size(Yf,1);

if previewFlag == 1
    Wp = data.Wp;
    Wf = data.Wf;
    wini = data.wini;
    wf = data.wf;
    pDist = size(Wp,1);
    fDist = size(Wf,1);
end

N = controlParams.N;
p = controlParams.p;
f = controlParams.f;
Q = controlParams.Q;
R = controlParams.R;
Nbar = N-p-f+1;

% Number of I/O channels
nInputs = pIn/p;
% nOutputs = pOut/p;

switch ivFlag
    case 0
        Z = eye(Nbar);
    case 1
        if previewFlag == 1
            Z = [Up;Yp;Uf;Wp;Wf];
        else
            Z = [Up;Yp;Uf];
        end
end

if method == 1 || method == 2  % QP or SDP
    % Decision variables
    g = ones(size(Z,1),1);
    uf = (pi/180)*ones(fIn,1);
    z0 = [g; uf];    

    % Equality constraints: Aeq*z = beq
    Aeq = [Up*Z' zeros(pIn,fIn);
        Yp*Z' zeros(pOut,fIn);
        Uf*Z' -eye(fIn)];

    beq = [uini;
        yini;
        zeros(fIn,1)];

    if previewFlag == 1
        Aeq = [Aeq;
            Wp*Z' zeros(pDist,fIn);
            Wf*Z' zeros(fDist,fIn)];

        beq = [beq;
            wini;
            wf];
    end

    % Inequality constraints: A*z <= b
    
    % Input rate conbstraints
    [A,b] = getInputRateConstr(nInputs,f,controlParams.duf);
    Aineq = [zeros(size(A,1),length(g)) A;
         zeros(size(A,1),length(g)) -A];
     
    bineq = [b(1:nInputs) + data.uini(end-nInputs+1:end); 
            b(nInputs+1:end);
            b(1:nInputs) - data.uini(end-nInputs+1:end);
            b(nInputs+1:end);];   

    % Input constraints
    A_lbu = -[zeros(fIn,length(g)) eye(fIn)];
    A_ubu = [zeros(fIn,length(g)) eye(fIn)];

    A_lby = -[Yf*Z' zeros(fOut,fIn)];
    A_uby = [Yf*Z' zeros(fOut,fIn)];

    Aineq = [Aineq;
             A_lbu;
             A_ubu;
             A_lby;
             A_uby];

    b_lbu = -kron(ones(f,1), controlParams.lbu);
    b_ubu = kron(ones(f,1), controlParams.ubu);

    b_lby = -kron(ones(f,1), controlParams.lby);
    b_uby = kron(ones(f,1), controlParams.uby);

    bineq = [bineq;
             b_lbu;
             b_ubu;
             b_lby;
             b_uby]; 

    if method == 1 % QP

    %%%%%%%%%%%%%%%%%%%% Constraints feasibility %%%%%%%%%%%%%%%%%%%%%%%%%%
    % [ineq_iis,eq_iis,~] = deletionfilter(Aineq,bineq,Aeq,beq,[],[]);    
    % find(ineq_iis); Infeasible_ineqs = numel(find(ineq_iis))
    % find(eq_iis); Infeasible_eqs = numel(find(eq_iis))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% active-set quadprog feasibility %%%%%%%%%%%%%%%%%%%%%%
    % H = 2* blkdiag(Z*Yf'*Q*Yf*Z',R);
    % HnSp = (null(Aeq)).'*H*null(Aeq);
    % posdef = min(real(eig(HnSp)))>=0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = length(z0);
    z = optimvar('z',n);    
    objFcn = (Yf*Z'*z(1:length(g))-rf)'*Q*(Yf*Z'*z(1:length(g))-rf) + ...
        z(length(g)+1:end)'*R*z(length(g)+1:end);
    qprob = optimproblem("Objective",objFcn);

    constr.eq = Aeq*z == beq;
    constr.ineq = Aineq*z <= bineq;   
    qprob.Constraints = constr;

    opts = optimoptions('quadprog','Algorithm','active-set');
    % opts = optimoptions('quadprog','Algorithm','interior-point-convex');
    initPoint = struct('z',z0);
    [sol,~,~,~] = solve(qprob,initPoint,'options',opts);
    z = sol.z;

    %%%%%%%%%%%%%%%%%%%%%%%%%%% quadprog %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % H1 = Z*Yf.'*Q.^(1/2);
    % H1 = H1*H1';
    % H = 2* blkdiag(H1,R);
    % % H = 2* blkdiag(Z*Yf'*Q*Yf*Z',R);
    % c = [-2*Z*Yf'*Q*rf;
    %     zeros(fIn,1)];  
    % options = optimoptions('quadprog','Algorithm', ...
    %     'interior-point-convex','MaxIterations',500);
    % z = quadprog(H,c,Aineq,bineq,Aeq,beq,[],[],[],options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    %%%%%%%%%%%%%%%%%%%%%%%% yalmip quadprog %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % x = sdpvar(size(z0,1),1);
    % Aeqx = Aeq*x;
    % Aineqx = Aineq*x;
    % Constraints = [Aeqx == beq, Aineqx <= bineq];
    % H1 = Z*Yf.'*Q.^(1/2);
    % H = 2* blkdiag(H1*H1',R);    
    % c = [-2*Z*Yf'*Q*rf;
    %     zeros(fIn,1)];  
    % 
    % Objective = (1/2)*x.'*H*x + c.'*x;
    % optimize(Constraints, Objective);
    % z = value(x);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     %%%%%%%%%%%%%%%%%%%%%%%%%%%% fmincon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % objFcn = @(z) (Yf*Z'*z(1:length(g))-rf)'*Q*(Yf*Z'*z(1:length(g))-rf) + ...
    %     z(length(g)+1:end)'*R*z(length(g)+1:end);
    % 
    % options = optimoptions('fmincon','Algorithm','interior-point','UseParallel',true);
    % z = fmincon(objFcn,z0,Aineq,bineq,Aeq,beq,lbu,ubu,[],options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif method == 2 % SDP
    H1 = Z*Yf.'*Q.^(1/2); % force H symmetric
    H = 2* blkdiag(H1*H1',R);   
    c = [-2*Z*Yf'*Q*rf;
        zeros(fIn,1)];  

    [L, D, P] = modchol_ldlt(H);   % Get H = M'*M
    R = P'*L*D^(1/2);
    M = R';

    % Collapse all constraints
    % A = [Aeq;
    %     -Aeq;
    %     Aineq];
    % 
    % b = [beq;
    %     -beq;
    %     bineq];   

    A = Aineq;
    b = bineq;  

    % Optimization variables
    t = sdpvar(1);
    x = sdpvar(size(H,1),1);

    % Objective function
    obj = t;

    % Constraints
    C = [eye(size(M,1)) M*x zeros(size(M,1),size(A,1));
        x.'*M.' t-c.'*x zeros(size(A,1),1).';
        zeros(size(A,1),size(M,1)) zeros(size(A,1),1) diag(b-A*x)];
    Constr = [C >=0, Aeq*x == beq]; 

    % Solve optimization
    % optimize(C >= 0, obj,sdpsettings('debug',1,'solver','sedumi'));
    % optimize(C >= 0, obj,sdpsettings('debug',1,'solver','sdpt3'));
    % optimize(Constr, obj,sdpsettings('debug',1,'solver','mosek'));
    optimize(Constr, obj,sdpsettings('debug',1,'solver','sdpt3'));
    z = value(x);    
    end        
elseif method ==3 % casadi + NLP
    g = casadi.SX.sym('g', size(Z,1), 1);
    uf = casadi.SX.sym('uf',fIn,1);
    z = [g;uf];

    if previewFlag == 1
        Aeq = [Up;Yp;Uf;Wp;Wf]*Z';
        beq = [uini;yini;uf;wini;wf];
    else
        Aeq = [Up;Yp;Uf]*Z';
        beq = [uini;yini;uf];
    end

    q = Aeq*g - beq;

    lbq = zeros(length(q),1);
    ubq = zeros(length(q),1);

    % Input rate constraints
    [A,b] = getInputRateConstr(nInputs,f,controlParams.duf);
    q = [q;
         A*uf];    
    
    lbq = [lbq; 
            -b(1) + data.uini(end);
            -b(2:end)];
    ubq = [ubq; 
            b(1) + data.uini(end);
            b(2:end)];

    % Bounds on input: lbu <= u <= ubu
    lbu = -inf(size(z));
    ubu = inf(size(z));

    lbu(length(g)+1:length(g)+length(uf)) = kron(ones(f,1), ...
        controlParams.lbu);
    ubu(length(g)+1:length(g)+length(uf)) = kron(ones(f,1), ...
        controlParams.ubu);

    cost = (Yf*Z'*g-rf)'*Q*(Yf*Z'*g-rf) + uf'*R*uf;

    prob = struct('f', cost, 'x', z, 'g', q);
    solver = casadi.nlpsol('solver', 'ipopt', prob);
    % solver = casadi.qpsol('solver','qpoases', prob);
    output  = solver('x0', zeros(size(z,1),1),'lbx', lbu, 'ubx', ubu, ...
        'lbg', lbq,'ubg', ubq);

    z = full(output.x);
else
    disp('Method doesn''t exist')
    return;
end

% g = z(1:size(Z,1));
uf = z(size(Z,1)+1:end);
uOpt = uf(1:nInputs);
% uOpt = uf(1:2*nInputs);

end

