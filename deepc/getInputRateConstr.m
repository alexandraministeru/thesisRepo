function [Aineq,bineq] = getInputRateConstr(nInputs,f,duf)
%getInputRateConstr(nInputs,f,duf) Function used to construct the input
%rate constraint matrices such that Aineq*uf <= bineq, with Aineq*uf = 
% [u1-u0 u2-u1 u3-u2 ... uf-uf-1]' and u0 being the previous input command.
%
% Input arguments
%----------------
% nInputs : scalar specifying the number of input channels. 
% f       : the length of the future window.
% duf     : column vector containing input rate constraints for each input 
%           channel [du_1 du_2 ... du_nInputs]'.
% Output arguments:
%------------------
% Aineq   : f*nInputs-by-f*nIputs matrix of the form 
%           [ I  0  0  0 ... 0  0
%            -I  I  0  0 ... 0  0
%             0 -I  I  0 ... 0  0
%             ...................
%             0  0  0  0 ...-I  I].
% bineq   : f*nInputs-by-1 matrix.
%==========================================================================

a = eye(f) + diag(-ones(1,f-1),-1);
Aineq = kron(a,eye(nInputs));
bineq = kron(ones(f,1),duf);

end