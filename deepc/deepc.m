function uOpt = deepc(data,uini,yini,N,p,f,rf,controlParams,method,ivFlag)
%deepc(Up,Uf,Yp,Yf,uini,yini,N,p,f,rf,Q,R,method) Function used to compute
%the DeePC optimal input.
%
% Input arguments
%----------------
% data          : structure containing the Hankel matrices Up,Yp,Uf,Yf
%                 (past inputs, past outputs, future inputs, and future
%                 outputs respectively)
% uini          : column vector containing past inputs.
% yini          : column vector containing past outputs.
% N             : length of data set used for DeePC.
% p             : past window size.
% f             : future window size.
% rf            : vector of reference values for the next f steps.
% controlParams : structure containing the DeePC output weight Q, DeePC
%                 input weight R, input lower bound lbu, and input upper
%                 bound ubu.
% method        : Optimization method (1=fmincon+sqp, 2=quadprog).
% ivFlag        : Flag that indicates whether or not instrumental variables
%                 are used (1 - IV used, 0 - IV not used)
%
% Output arguments:
%------------------
% uOpt          : Optimal control input.
%==========================================================================
Nbar = N-p-f+1;
Up = data.Up;
Yp = data.Yp;
Uf = data.Uf;
Yf = data.Yf;
Q = controlParams.Q;
R = controlParams.R;

nInputs = size(Up,1)/p;
nOutputs = size(Yp,1)/p;

if ivFlag == 0 % no IV
    % Decision variables
    g = ones(Nbar,1);
    uf = ones(size(Uf,1),1);
    r = zeros(f*nOutputs,1);
    z = [g; uf; r];

    % Bounds on input: lbu <= u <= ubu
    lbu = -inf(size(z));
    ubu = inf(size(z));

    lbu(length(g)+1:length(g)+length(uf)) = kron(ones(f,1), ...
        controlParams.lbu);
    ubu(length(g)+1:length(g)+length(uf)) = kron(ones(f,1), ...
        controlParams.ubu);

    % Equality constraints: Aeq*z = beq
    Aeq = [Up zeros(size(Up,1),size(Uf,1)) zeros(size(Up,1),size(Yf,1));
        Yp zeros(size(Up,1),size(Uf,1)) zeros(size(Up,1),size(Yf,1));
        Uf -eye(size(Uf,1)) zeros(size(Uf,1),size(Yf,1));
        zeros(size(Yf,1),Nbar) zeros(size(Yf,1),size(Uf,1)) eye(size(Yf,1))];

    beq = [uini;
        yini;
        zeros(size(Uf,1),1);
        rf];

    if method == 1 % fmincon
        % yf = Yf*g
        % min(uf,g) (yf-rf)'*Q*(yf-rf) + uf'*R*uf <=> min(uf,g) (g'*Yf'-rf')*Q*(Yf*g-rf) + uf'*R*uf
        objFcn = @(z) (z' * [Yf'*Q*Yf zeros(Nbar,size(Uf,1)) -Yf'*Q; ...
            zeros(size(Uf,1),Nbar) R zeros(size(Uf,1),size(Yf,1)); ...
            -Q*Yf zeros(size(Yf,1),size(Uf,1)) Q] * z);
        options = optimoptions('fmincon','Algorithm','sqp');
        z = fmincon(objFcn,z,[],[],Aeq,beq,lbu,ubu,[],options);

    elseif method == 2 % quadprog
        H = 2*[Yf'*Q*Yf zeros(Nbar,size(Uf,1)) -Yf'*Q; ...
            zeros(size(Uf,1),Nbar) R zeros(size(Uf,1),size(Yf,1)); ...
            -Q*Yf zeros(size(Yf,1),size(Uf,1)) Q];
        options = optimoptions('quadprog', 'MaxIterations',200);
        z = quadprog(H,[],[],[],Aeq,beq,lbu,ubu,[],options);

    elseif method == 3 % casadi
        g = casadi.SX.sym('g', size([Up;Yp;Uf],2), 1);
        uf = casadi.SX.sym('uf',size(Uf,1),1);
        z = [g;uf];

        Aeq = [Up;Yp;Uf];
        beq = [uini;yini;uf];

        q = Aeq*g - beq;
        cost = (Yf*g-rf)'*Q*(Yf*g-rf) + uf'*R*uf;

        lbq = zeros(length(q),1);
        ubq = zeros(length(q),1);

        lbu = lbu(1:size(z,1));
        ubu = ubu(1:size(z,1));

        prob = struct('f', cost, 'x', z, 'g', q);
        solver = casadi.nlpsol('solver', 'ipopt', prob);
        output  = solver('x0', zeros(size(z,1),1),'lbx', lbu, 'ubx', ubu, ...
            'lbg', lbq,'ubg', ubq);

        z = full(output.x);

    elseif method == 4 % casadi QP
        H = 2*casadi.DM([Yf'*Q*Yf zeros(Nbar,size(Uf,1)) -Yf'*Q; ...
            zeros(size(Uf,1),Nbar) R zeros(size(Uf,1),size(Yf,1)); ...
            -Q*Yf zeros(size(Yf,1),size(Uf,1)) Q]);

        g = casadi.DM.zeros(size(H,1),1);

        Aeq = casadi.DM([Up zeros(size(Up,1),size(Uf,1)) zeros(size(Up,1),size(Yf,1));
            Yp zeros(size(Up,1),size(Uf,1)) zeros(size(Up,1),size(Yf,1));
            Uf -eye(size(Uf,1)) zeros(size(Uf,1),size(Yf,1));
            zeros(size(Yf,1),Nbar) zeros(size(Yf,1),size(Uf,1)) eye(size(Yf,1))]);

        beq = casadi.DM([uini;
            yini;
            zeros(size(Uf,1),1);
            rf]);

        prob = struct;
        prob.h = H.sparsity();
        prob.a = Aeq.sparsity();
        solver = casadi.conic('solver','qpoases',prob);

        output = solver('h', H, 'g', g, 'a', Aeq, 'lba', beq, 'uba', beq,'lbx',lbu,'ubx',ubu);
        z = full(output.x);
    end
    % g = z(1:Nbar);
    uf = z(length(g)+1:end);

else % with IV
    % Instrumental variables
    Z = [Up;Yp;Uf];

    % Decision variables
    g_hat = zeros(size(Z,1),1);
    uf = zeros(size(Uf,1),1);
    r = zeros(f*nOutputs,1);
    z = [g_hat; uf; r];

    % Bounds on input: lbu <= u <= ubu
    lbu = -inf(size(z));
    ubu = inf(size(z));

    lbu(length(g_hat)+1:length(g_hat)+length(uf)) = kron(ones(f,1), ...
        controlParams.lbu);
    ubu(length(g_hat)+1:length(g_hat)+length(uf)) = kron(ones(f,1), ...
        controlParams.ubu);

    % Equality constraints: Aeq*z = beq
    Aeq = [Up*Z' zeros(size(Up,1),size(Uf,1)) zeros(size(Up,1),size(Yf,1));
        Yp*Z' zeros(size(Yp,1),size(Uf,1)) zeros(size(Yp,1),size(Yf,1));
        Uf*Z' -eye(size(Uf,1)) zeros(size(Uf,1),size(Yf,1));
        zeros(size(Yf,1),size(Z,1)) zeros(size(Yf,1),size(Uf,1)) eye(size(Yf,1))];

    beq = [uini;
        yini;
        zeros(size(Uf,1),1);
        rf];

    if method == 1 % fmincon
        % yf = Yf*g
        % min(uf,g) (yf-rf)'*Q*(yf-rf) + uf'*R*uf <=> min(uf,g) (g'*Yf'-rf')*Q*(Yf*g-rf) + uf'*R*uf
        objFcn = @(z) (z' * [Z*Yf'*Q*Yf*Z' zeros(size(Z,1),size(Uf,1)) -Z*Yf'*Q; ...
            zeros(size(Uf,1),size(Z,1)) R zeros(size(Uf,1),size(Yf,1)); ...
            -Q*Yf*Z' zeros(size(Yf,1),size(Uf,1)) Q] * z);
        options = optimoptions('fmincon','Algorithm','sqp');
        z = fmincon(objFcn,z,[],[],Aeq,beq,lbu,ubu,[],options);

    elseif method == 2 % quadprog
        H = 2*[Z*Yf'*Q*Yf*Z' zeros(size(Z,1),size(Uf,1)) -Z*Yf'*Q; ...
            zeros(size(Uf,1),size(Z,1)) R zeros(size(Uf,1),size(Yf,1)); ...
            -Q*Yf*Z' zeros(size(Yf,1),size(Uf,1)) Q];
        options = optimoptions('quadprog', 'MaxIterations',200);
        z = quadprog(H,[],[],[],Aeq,beq,lbu,ubu,[],options);

    elseif method == 3 % casadi
        g_hat = casadi.SX.sym('g_hat', size(Z,1), 1);
        uf = casadi.SX.sym('uf',size(Uf,1),1);
        z = [g_hat;uf];

        Aeq = [Up;Yp;Uf]*Z';
        beq = [uini;yini;uf];

        q = Aeq*g_hat - beq;
        cost = (Yf*Z'*g_hat-rf)'*Q*(Yf*Z'*g_hat-rf) + uf'*R*uf;

        lbq = zeros(length(q),1);
        ubq = zeros(length(q),1);

        lbu = lbu(1:size(z,1));
        ubu = ubu(1:size(z,1));

        prob = struct('f', cost, 'x', z, 'g', q);
        solver = casadi.nlpsol('solver', 'ipopt', prob);
        % solver = casadi.qpsol('solver','qpoases', prob);
        output  = solver('x0', zeros(size(z,1),1),'lbx', lbu, 'ubx', ubu, ...
            'lbg', lbq,'ubg', ubq);

        z = full(output.x);

    elseif method == 4 % casadi QP
        H = 2*casadi.DM([Z*Yf'*Q*Yf*Z' zeros(size(Z,1),f) -Z*Yf'*Q; ...
            zeros(f,size(Z,1)) R zeros(f,f); ...
            -Q*Yf*Z' zeros(f,f) Q]);

        g = casadi.DM.zeros(size(H,1),1);

        Aeq = casadi.DM([Up*Z' zeros(size(Up,1),size(Uf,1)) zeros(size(Up,1),size(Yf,1));
            Yp*Z' zeros(size(Yp,1),size(Uf,1)) zeros(size(Yp,1),size(Yf,1));
            Uf*Z' -eye(size(Uf,1)) zeros(size(Uf,1),size(Yf,1));
            zeros(size(Yf,1),size(Z,1)) zeros(size(Yf,1),size(Uf,1)) eye(size(Yf,1))]);

        beq = casadi.DM([uini;
            yini;
            zeros(size(Uf,1),1);
            rf]);

        prob = struct;
        prob.h = H.sparsity();
        prob.a = Aeq.sparsity();
        solver = casadi.conic('solver','qpoases',prob);

        output = solver('h', H, 'g', g, 'a', Aeq, 'lba', beq, 'uba', beq,'lbx',lbu,'ubx',ubu);
        z = full(output.x);
    end
    % g_hat = z(1:size(Z,1));
    uf = z(size(Z,1)+1:end);
end

uOpt = uf(1:nInputs);

end