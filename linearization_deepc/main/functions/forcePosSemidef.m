function H_pos_semidef = forcePosSemidef(H)
[V,d] = eig(H,'vector');
tol = length(d)*eps(max(d));
% eigs((real(d)<0)&(real(d)>-tol)) = 0;

if(find((real(d)<0)&(real(d)<-tol)))
    disp('Eigenvalues with real part lower than the tolerance were forced to zero');
end
d(d<0) = 0;
H_pos_semidef = real(V*diag(d)/V);
end