function uOpt = deepc(data,uini,yini,N,p,f,rf,Q,R,method,ivFlag)
%deepc(Up,Uf,Yp,Yf,uini,yini,N,p,f,rf,Q,R,method) Function used to compute
%the DeePC optimal input.
%
% Input arguments
%----------------
% data  : structure containing the Hankel matrices Up,Yp,Uf,Yf (past
%         inputs, past outputs, future inputs, and future outputs
%         respectively)
% uini  : column vector containing past inputs.
% yini  : column vector containing past outputs.
% N     : length of data set used for DeePC.
% p     : past window size.
% f     : future window size.
% rf    : vector of reference values for the next f steps.
% Q     : DeePC output weight.
% R     : DeePC input weight.
% method: Optimization method (1=fmincon+sqp, 2=quadprog).
% ivFlag: Flag that indicates whether or not instrumental variables are
%         used (1 - IV used, 0 - IV not used)
%
% Output arguments:
%------------------
% uOpt  : Optimal control input.
%==========================================================================
Nbar = N-p-f+1;
% Decision variables
g = ones(Nbar,1);
uf = ones(f,1);
r = zeros(f,1);
z = [g; uf; r];

Up = data.Up;
Yp = data.Yp;
Uf = data.Uf;
Yf = data.Yf;

if ivFlag == 0 % no IV
    % Decision variables
    g = ones(Nbar,1);
    uf = ones(f,1);
    r = zeros(f,1);
    z = [g; uf; r];

    % Equality constraints: Aeq*z = beq
    Aeq = [Up zeros(p,f) zeros(p,f);
        Yp zeros(p,f) zeros(p,f);
        Uf -eye(f,f) zeros(f,f);
        zeros(f,Nbar) zeros(f,f) eye(f,f)];
    beq = [uini;
        yini;
        zeros(f,1);
        rf];

    if method == 1 % fmincon
        % yf = Yf*g
        % min(uf,g) (yf-rf)'*Q*(yf-rf) + uf'*R*uf <=> min(uf,g) (g'*Yf'-rf')*Q*(Yf*g-rf) + uf'*R*uf
        objFcn = @(z) (z' * [Yf'*Q*Yf zeros(Nbar,f) -Yf'*Q; ...
            zeros(f,Nbar) R zeros(f,f); ...
            -Q*Yf zeros(f,f) Q] * z);
        options = optimoptions('fmincon','Algorithm','sqp');
        z = fmincon(objFcn,z,[],[],Aeq,beq,[],[],[],options);

    elseif method == 2 % quadprog
        H = 2*[Yf'*Q*Yf zeros(Nbar,f) -Yf'*Q; ...
            zeros(f,Nbar) R zeros(f,f); ...
            -Q*Yf zeros(f,f) Q];
        z = quadprog(H,[],[],[],Aeq,beq);

    elseif method == 3 % casadi
        g = casadi.SX.sym('g', size([Up;Yp;Uf],2), 1);
        uf = casadi.SX.sym('uf',size(Uf,1),1);
        z = [g;uf];

        Aeq = [Up;Yp;Uf];
        beq = [uini;yini;uf];

        q = Aeq*g - beq;
        cost = (Yf*g-rf)'*Q*(Yf*g-rf) + uf'*R*uf;

        lbq     = zeros(length(q),1);
        ubq     = zeros(length(q),1);

        prob = struct('f', cost, 'x', z, 'g', q);
        solver = nlpsol('solver', 'ipopt', prob);
        output  = solver('x0', zeros(size(z,1),1),'lbg', lbq,'ubg', ubq);

        z = full(output.x);
    end
    g = z(1:Nbar);
    uf = z(Nbar+1:end);

else % with IV
     % Instrumental variables
    Z = [Up;Yp;Uf];

    % Decision variables
    g_hat = ones(size(Z,1),1);
    uf = ones(f,1);
    r = zeros(f,1);
    z = [g_hat; uf; r];
   
    % Equality constraints: Aeq*z = beq
    Aeq = [Up*Z' zeros(p,f) zeros(p,f);
        Yp*Z' zeros(p,f) zeros(p,f);
        Uf*Z' -eye(f,f) zeros(f,f);
        zeros(f,size(Z,1)) zeros(f,f) eye(f,f)];

    beq = [uini;
        yini;
        zeros(f,1);
        rf];

    if method == 1 % fmincon
        % yf = Yf*g
        % min(uf,g) (yf-rf)'*Q*(yf-rf) + uf'*R*uf <=> min(uf,g) (g'*Yf'-rf')*Q*(Yf*g-rf) + uf'*R*uf
        objFcn = @(z) (z' * [Z*Yf'*Q*Yf*Z' zeros(size(Z,1),f) -Z*Yf'*Q; ...
            zeros(f,size(Z,1)) R zeros(f,f); ...
            -Q*Yf*Z' zeros(f,f) Q] * z);
        options = optimoptions('fmincon','Algorithm','sqp');
        z = fmincon(objFcn,z,[],[],Aeq,beq,[],[],[],options);
    elseif method == 2 % quadprog
        H = 2*[Z*Yf'*Q*Yf*Z' zeros(size(Z,1),f) -Z*Yf'*Q; ...
            zeros(f,size(Z,1)) R zeros(f,f); ...
            -Q*Yf*Z' zeros(f,f) Q];
        z = quadprog(H,[],[],[],Aeq,beq);
    elseif method == 3 % casadi
        %TBD
    end
    % g_hat = z(1:size(Z,1));
    uf = z(size(Z,1)+1:end);
end

uOpt = uf(1); %HARDCODED INPUT SIZE

end