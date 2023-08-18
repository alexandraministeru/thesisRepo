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
% method        : Optimization method (1=fmincon+sqp, 2=quadprog, 3=casadi).
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

% Input/output samples within past/future windows
pIn = size(Up,1);
pOut = size(Yp,1);
fIn = size(Uf,1);
fOut = size(Yf,1);

nInputs = size(Up,1)/p;
% nOutputs = size(Yp,1)/p;

switch ivFlag
    case 0
        Z = eye(Nbar);
    case 1
        Z = [Up;Yp;Uf];
end

if (method == 1 || method == 2) % fmincon or quadprog
    % Decision variables
    g = ones(size(Z,1),1);
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

    if method == 1 % fmincon
        % yf = Yf*g
        % min(uf,g) (yf-rf)'*Q*(yf-rf) + uf'*R*uf <=> min(uf,g) (g'*Yf'-rf')*Q*(Yf*g-rf) + uf'*R*uf
        objFcn = @(z) (z' * [Z*Yf'*Q*Yf*Z' zeros(size(Z,1),fIn) -Z*Yf'*Q; ...
            zeros(fIn,size(Z,1)) R zeros(fIn,fOut); ...
            -Q*Yf*Z' zeros(fOut,fIn) Q] * z);
        options = optimoptions('fmincon','Algorithm','sqp');
        z = fmincon(objFcn,z,[],[],Aeq,beq,lbu,ubu,[],options);

    elseif method == 2 % quadprog
        H = 2*[Z*Yf'*Q*Yf*Z' zeros(size(Z,1),fIn) -Z*Yf'*Q; ...
            zeros(fIn,size(Z,1)) R zeros(fIn,fOut); ...
            -Q*Yf*Z' zeros(fOut,fIn) Q];
        options = optimoptions('quadprog', 'MaxIterations',200);
        z = quadprog(H,[],[],[],Aeq,beq,lbu,ubu,[],options);
    end

elseif method == 3 % casadi
    g = casadi.SX.sym('g', size(Z,1), 1);
    uf = casadi.SX.sym('uf',fIn,1);
    z = [g;uf];

    Aeq = [Up;Yp;Uf]*Z';
    beq = [uini;yini;uf];

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