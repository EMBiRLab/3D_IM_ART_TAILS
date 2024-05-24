function [f, gradf] = objective_fun(z, nsteps, nstates, dt, xref, PRY_idx_in_q, varargin)

    if nargin == 7
        coeff = varargin{1};
    else
        coeff = 1;
    end
    
    ind_mtx = sparse(length(PRY_idx_in_q), nstates/2);
    ind_mtx(:, PRY_idx_in_q) = eye(length(PRY_idx_in_q));

    %%
    X = reshape(z(1:(nsteps+1)*nstates, 1), nstates, []);
    
    q = X(1:nstates/2, :);
    qd = X(nstates/2+1:end, :);
    q_mid = (q(:, 1:nsteps) + q(:, 2:nsteps+1))/2 + (qd(:, 1:nsteps) - qd(:, 2:nsteps+1))*dt/8;
    
    diff4f1 = angdiff(reshape(xref(:, 1:2:2*nsteps-1), [], 1), reshape(ind_mtx*q(:, 1:nsteps), [], 1));
    f1 = diff4f1'*diff4f1;

    diff4f2 = angdiff(reshape(xref(:, 2:2:2*nsteps), [], 1), reshape(ind_mtx*q_mid, [], 1));
    f2 = 4*(diff4f2'*diff4f2);

    diff4f3 = angdiff(reshape(xref(:, 3:2:2*nsteps+1), [], 1), reshape(ind_mtx*q(:, 2:nsteps+1), [], 1));
    f3 = diff4f3'*diff4f3;

    f = (f1 + f2 + f3)*dt/6;
    
    %%
    gradf1 = sparse(1, length(z));

    idx1_mtx = (1:nstates:(nsteps-1)*nstates+1)' + linspace(0, nstates/2-1, nstates/2);
    idx1_vec = reshape(idx1_mtx', 1, []);

    gradf1(1, idx1_vec) = reshape((q(:, 1:nsteps)'*(2*(ind_mtx'*ind_mtx)) - 2*xref(:, 1:2:2*nsteps-1)'*ind_mtx)', 1, []);

    %%
    gradf2_comp1 = sparse(1, length(z));
    gradf2_comp2 = sparse(1, length(z));
    
    gradf2_comp1(1, 1:nsteps*nstates) = reshape((4*(q_mid'*(2*(ind_mtx'*ind_mtx)) - 2*xref(:, 2:2:2*nsteps)'*ind_mtx) * ...
                                                [eye(nstates/2)/2, eye(nstates/2)*dt/8])', 1, []);

    gradf2_comp2(1, nstates+1:(nsteps+1)*nstates) = reshape((4*(q_mid'*(2*(ind_mtx'*ind_mtx)) - 2*xref(:, 2:2:2*nsteps)'*ind_mtx) * ...
                                                [eye(nstates/2)/2, -eye(nstates/2)*dt/8])', 1, []);
    gradf2 = gradf2_comp1 + gradf2_comp2;

    %%
    gradf3 = sparse(1, length(z));

    idx3_mtx = (nstates+1:nstates:nsteps*nstates+1)' + linspace(0, nstates/2-1, nstates/2);
    idx3_vec = reshape(idx3_mtx', 1, []);

    gradf3(1, idx3_vec) = reshape((q(:, 2:nsteps+1)'*(2*(ind_mtx'*ind_mtx)) - 2*xref(:, 3:2:2*nsteps+1)'*ind_mtx)', 1, []);

    %%
    gradf = (gradf1 + gradf2 + gradf3)*dt/6;
    
    %%
    f = coeff*f;
    gradf = coeff*gradf;

end

