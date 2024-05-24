function [Aeq, beq] = linear_eq_con(nsteps, nstates, ninputs, x0, link_num, total_len)
    
    ndec = (nsteps + 1)*nstates + (nsteps*2 + 1)*ninputs + link_num;
    
    % initial condition equality constraint
    A_x0 = sparse(nstates, ndec);
    b_x0 = sparse(nstates, 1);

    A_x0(1 : nstates, 1 : nstates) = eye(nstates);
    b_x0(1 : nstates, 1) = x0;
    
    % total length equality constraint
    A_totallen = sparse(1, ndec);
    
    A_totallen(1, ndec-link_num+1 : end) = ones(1, link_num);
    b_totallen = total_len;
    
    % stack all the constraints
    Aeq = [A_x0; A_totallen];
    beq = [b_x0; b_totallen];

end
