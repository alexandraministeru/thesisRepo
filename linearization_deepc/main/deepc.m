function uOpt = deepc(data,rf,controlParams,method,ivFlag,previewFlag)
% deepc(data,rf,controlParams,method,ivFlag,previewFlag) Function used to 
% compute the DeePC optimal input.
%
% Input arguments
%----------------
% data          : structure containing the Hankel matrices
%                 Up,Yp,Uf,Yf,Wp,Wf (past inputs, past outputs, future
%                 inputs, future outputs, past disturbance, and future
%                 disturbance respectively), uini,yini,wini,wf (past
%                 input, output and distubance trajectory and disturbance
%                 preview respectively). Wp,Wf,wini,wf are not used if
%                 previewFlag is set to 0.
% rf            : vector of reference values for the next f steps.
% controlParams : structure containing the DeePC output weight Q, DeePC
%                 input weight R, input lower bound lbu, and input upper
%                 bound ubu, the length of the data set used for DeePC N, 
%                 past window size p, and future window size f.
% method        : Optimization method (1 - fmincon+sqp, 2 - quadprog, 
%                 3 - casadi+nlp, 4 - casadi+qp).
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

if (method == 1 || method == 2) % fmincon or quadprog
    % Decision variables
    g = zeros(size(Z,1),1);
    uf = ones(fIn,1);
    r = zeros(fOut,1);
    z = [g; uf; r];

    % Bounds on input: lbu <= u <= ubu
    lbu = -inf(size(z));
    ubu = inf(size(z));

    lbu(length(g)+1:length(g)+length(uf)) = kron(ones(f,1), ...
        controlParams.lbu);
    ubu(length(g)+1:length(g)+length(uf)) = kron(ones(f,1), ...
        controlParams.ubu);

    % Equality constraints: Aeq*z = beq

    Aeq = [Up*Z' zeros(pIn,fIn) zeros(pIn,fOut);
        Yp*Z' zeros(pOut,fIn) zeros(pOut,fOut);
        Uf*Z' -eye(fIn) zeros(fIn,fOut);
        zeros(fOut,size(Z,1)) zeros(fOut,fIn) eye(fOut)];

    beq = [uini;
        yini;
        zeros(fIn,1);
        rf];

    if previewFlag == 1
        Aeq = [Aeq;
            Wp*Z' zeros(pDist,fIn) zeros(pDist,fOut);
            Wf*Z' zeros(fDist,fIn) zeros(fDist,fOut)];

        beq = [beq;
            wini;
            wf];
    end

    if method == 1 % fmincon
        % yf = Yf*g
        % min(uf,g) (yf-rf)'*Q*(yf-rf) + uf'*R*uf <=>
        % min(uf,g) (g'*Yf'-rf')*Q*(Yf*g-rf) + uf'*R*uf
        % s.t Aeq*z = beq, lbu <= uf <= ubu

        objFcn = @(z) (z' * [Z*Yf'*Q*Yf*Z' zeros(size(Z,1),fIn) -Z*Yf'*Q; ...
            zeros(fIn,size(Z,1)) R zeros(fIn,fOut); ...
            -Q*Yf*Z' zeros(fOut,fIn) Q] * z);
        options = optimoptions('fmincon','Algorithm','sqp');
        options = optimoptions(options,'UseParallel',true);

        z = fmincon(objFcn,z,[],[],Aeq,beq,lbu,ubu,[],options);

    else % quadprog
        % min(uf,g) (1/2)z'Hz
        % s.t Aeq*z = beq,  lbu <= uf <= ubu

        H = 2*[Z*Yf'*Q*Yf*Z' zeros(size(Z,1),fIn) -Z*Yf'*Q; ...
            zeros(fIn,size(Z,1)) R zeros(fIn,fOut); ...
            -Q*Yf*Z' zeros(fOut,fIn) Q];
        options = optimoptions('quadprog', 'MaxIterations',200);

        z = quadprog(H,[],[],[],Aeq,beq,lbu,ubu,[],options);
    end

elseif method == 3 % casadi + NLP
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
    cost = (Yf*Z'*g-rf)'*Q*(Yf*Z'*g-rf) + uf'*R*uf;

    lbq = zeros(length(q),1);
    ubq = zeros(length(q),1);

    % Bounds on input: lbu <= u <= ubu
    lbu = -inf(size(z));
    ubu = inf(size(z));

    lbu(length(g)+1:length(g)+length(uf)) = kron(ones(f,1), ...
        controlParams.lbu);
    ubu(length(g)+1:length(g)+length(uf)) = kron(ones(f,1), ...
        controlParams.ubu);

    prob = struct('f', cost, 'x', z, 'g', q);
    solver = casadi.nlpsol('solver', 'ipopt', prob);
    % solver = casadi.qpsol('solver','qpoases', prob);
    output  = solver('x0', zeros(size(z,1),1),'lbx', lbu, 'ubx', ubu, ...
        'lbg', lbq,'ubg', ubq);

    z = full(output.x);

elseif method == 4 % casadi + QP
    H = 2*casadi.DM([Z*Yf'*Q*Yf*Z' zeros(size(Z,1),f) -Z*Yf'*Q; ...
        zeros(f,size(Z,1)) R zeros(f,f); ...
        -Q*Yf*Z' zeros(f,f) Q]);

    g = casadi.DM.zeros(size(H,1),1);

    if previewFlag == 1
        Aeq = casadi.DM([Up*Z' zeros(pIn,fIn) zeros(pIn,fOut);
        Yp*Z' zeros(pOut,fIn) zeros(pOut,fOut);
        Uf*Z' -eye(fIn) zeros(fIn,fOut);
        zeros(fOut,size(Z,1)) zeros(fOut,fIn) eye(fOut);
        Wp*Z' zeros(pDist,fIn) zeros(pDist,fOut);
        Wf*Z' zeros(fDist,fIn) zeros(fDist,fOut)]);

        beq = casadi.DM([uini;
            yini;
            zeros(fIn,1);
            rf;
            wini;
            wf;]);
    else
        Aeq = casadi.DM([Up*Z' zeros(pIn,fIn) zeros(pIn,fOut);
            Yp*Z' zeros(pOut,fIn) zeros(pOut,fOut);
            Uf*Z' -eye(fIn) zeros(fIn,fOut);
            zeros(fOut,size(Z,1)) zeros(fOut,fIn) eye(fOut)]);

        beq = casadi.DM([uini;
            yini;
            zeros(fIn,1);
            rf]);
    end

    % Bounds on input: lbu <= u <= ubu
    lbu = -inf(size(Z,1)+fIn+fOut,1);
    ubu = inf(size(Z,1)+fIn+fOut,1);

    lbu(length(g)+1:length(g)+fIn) = kron(ones(f,1), ...
        controlParams.lbu);
    ubu(length(g)+1:length(g)+fIn) = kron(ones(f,1), ...
        controlParams.ubu);

    prob = struct;
    prob.h = H.sparsity();
    prob.a = Aeq.sparsity();
    solver = casadi.conic('solver','qpoases',prob);

    output = solver('h', H, 'g', g, 'a', Aeq, 'lba', beq, 'uba', beq, ...
        'lbx',lbu,'ubx',ubu);

    z = full(output.x);
else
    disp('Method doesn''t exist')
    return;
end

% g = z(1:size(Z,1));
uf = z(size(Z,1)+1:end);

uOpt = uf(1:nInputs);

end

