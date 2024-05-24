function [Aeq, beq] = linear_eq_con(nsteps, nstates, ninputs, x0)
    
    ndec = (nsteps + 1)*nstates + (nsteps*2 + 1)*ninputs;
    
    A_x0 = sparse(nstates, ndec);
    b_x0 = sparse(nstates, 1);

    A_x0(1 : nstates, 1 : nstates) = eye(nstates);
    b_x0(1 : nstates, 1) = x0;

    Aeq = A_x0;
    beq = b_x0;

end
