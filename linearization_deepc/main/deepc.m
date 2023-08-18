function uOpt = deepc(data,rf,controlParams,method,ivFlag,previewFlag)
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
% fOut = size(Yf,1);

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

if method == 1 % quadprog
    % Decision variables
    g = ones(size(Z,1),1);
    uf = (pi/180)*ones(fIn,1);
    z0 = [g; uf];

    % Bounds on input: lbu <= u <= ubu
    lbu = -inf(size(z0));
    ubu = inf(size(z0));

    lbu(length(g)+1:length(g)+length(uf)) = kron(ones(f,1), ...
        controlParams.lbu);
    ubu(length(g)+1:length(g)+length(uf)) = kron(ones(f,1), ...
        controlParams.ubu);

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
    
    n = length(z0);
    z = optimvar('z',n,'LowerBound',lbu,'UpperBound',ubu);

    constr = Aeq*z == beq;

    objFcn = (Yf*Z'*z(1:length(g))-rf)'*Q*(Yf*Z'*z(1:length(g))-rf) + ...
        z(length(g)+1:end)'*R*z(length(g)+1:end);
    qprob = optimproblem("Objective",objFcn);
    qprob.Constraints = constr;
    opts = optimoptions('quadprog','Algorithm','interior-point-convex');
    z00 = struct('z',z0);

    [sol,~,~,~] = solve(qprob,z00,'options',opts);
    z = sol.z;
elseif method == 2 % casadi + NLP
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
else
    disp('Method doesn''t exist')
    return;
end

% g = z(1:size(Z,1));
uf = z(size(Z,1)+1:end);

uOpt = uf(1:nInputs);

end

