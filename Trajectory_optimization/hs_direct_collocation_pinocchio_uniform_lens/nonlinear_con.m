function [c, ceq, gradc, gradceq] = nonlinear_con(z, nsteps, nstates, ninputs, torso_dof, dt, power_limit, state_range,...
                                                  urdf_filename_drd, urdf_filepath_pino, varargin)

    if size(z, 2) > size(z, 1)
        z = z';
    end

    if nargin == 11
        coeff = varargin{1};
    else
        coeff.dyn = 1;
        coeff.pow = 1;
        coeff.coll = 1;
        coeff.xmid = 1;
    end
    
    ind_mtx = sparse(ninputs, nstates/2);
    ind_mtx(:, torso_dof+1:end) = eye(ninputs);

    xmid_ub = state_range(torso_dof+1:nstates/2, 2);
    xmid_lb = state_range(torso_dof+1:nstates/2, 1);

    xdmid_ub = state_range(nstates/2+torso_dof+1:nstates, 2);
    xdmid_lb = state_range(nstates/2+torso_dof+1:nstates, 1);

    %%
    dyn_ceq = zeros(nsteps*nstates, 1);
    dyn_gradceq = sparse(nsteps*nstates, length(z));

    pow_c = zeros(2*nsteps+1, 1);
    pow_gradc = sparse(2*nsteps+1, length(z));
    
    switch urdf_filename_drd
        case {'three-link_2JDOF_tail_robot_drd_3DOF_float.urdf'}
            coll_c = zeros(3*nsteps, 1);
            coll_gradc = sparse(3*nsteps, length(z));

        case {'four-link_2JDOF_tail_robot_drd_3DOF_float.urdf'}
            coll_c = zeros(6*nsteps, 1);
            coll_gradc = sparse(6*nsteps, length(z));
        
        case {'five-link_2JDOF_tail_robot_drd_3DOF_float.urdf'}
            coll_c = zeros(9*nsteps, 1);
            coll_gradc = sparse(9*nsteps, length(z));
        
        case {'six-link_2JDOF_tail_robot_drd_3DOF_float.urdf'}
            coll_c = zeros(12*nsteps, 1);
            coll_gradc = sparse(12*nsteps, length(z));
        
        case {'one-link_2JDOF_tail_robot_drd_3DOF_float.urdf', 'two-link_2JDOF_tail_robot_drd_3DOF_float.urdf'}
            coll_c = [];
            coll_gradc = [];
            
        otherwise
            error('Undefined urdf')
    end
    
    xmid_c = zeros(4*nsteps*ninputs, 1);
    xmid_gradc = sparse(4*nsteps*ninputs, length(z));

    %%
    X = reshape(z(1:(nsteps+1)*nstates, 1), nstates, []);
    
    Uact = reshape(z((nsteps+1)*nstates+1:end, 1), ninputs, []);
    Uact_mid = Uact(:, 2:2:2*nsteps);
    
    Uall = [zeros(torso_dof, 2*nsteps+1); Uact];
    Uall_mid = [zeros(torso_dof, nsteps); Uact_mid];

    q = X(1:nstates/2, :);
    qd = X(nstates/2+1:nstates, :);
    
    [qdd, p_qdd_p_q, p_qdd_p_qd, p_qdd_p_Uall] = differentiate_FDab_mex(urdf_filepath_pino, q, qd, Uall(:, 1:2:2*nsteps+1));
    
    q_mid = (q(:, 1:nsteps) + q(:, 2:nsteps+1))/2 + (qd(:, 1:nsteps) - qd(:, 2:nsteps+1))*dt/8;
    qd_mid = (qd(:, 1:nsteps) + qd(:, 2:nsteps+1))/2 + (qdd(:, 1:nsteps) - qdd(:, 2:nsteps+1))*dt/8;
    
    [qdd_mid, p_qdd_mid_p_q_mid, p_qdd_mid_p_qd_mid, p_qdd_mid_p_Uall_mid] = differentiate_FDab_mex(urdf_filepath_pino, ...
                                                                                                    q_mid, qd_mid, Uall_mid); 
    
    pos_ceq = q(:, 2:nsteps+1) - q(:, 1:nsteps) - (qd(:, 1:nsteps) + 4*qd_mid + qd(:, 2:nsteps+1))*dt/6;
    vel_ceq = qd(:, 2:nsteps+1) - qd(:, 1:nsteps) - (qdd(:, 1:nsteps) + 4*qdd_mid + qdd(:, 2:nsteps+1))*dt/6;

    xmid_c = reshape([ind_mtx*q_mid - xmid_ub; 
                      ind_mtx*qd_mid - xdmid_ub;
                      -ind_mtx*q_mid + xmid_lb;
                      -ind_mtx*qd_mid + xdmid_lb], [], 1);

    for i = 1 : nsteps
        dyn_ceq(nstates*(i-1)+1 : nstates*i, 1) = [pos_ceq(:, i); vel_ceq(:, i)];

        p_q_mid_p_q = eye(nstates/2)/2;
        p_q_mid_p_qd = eye(nstates/2)*dt/8;
        p_q_mid_p_qnext = eye(nstates/2)/2;
        p_q_mid_p_qdnext = -eye(nstates/2)*dt/8;

        p_qd_mid_p_q = (dt/8)*p_qdd_p_q(:, :, i);
        p_qd_mid_p_qd = eye(nstates/2)/2 + (dt/8)*p_qdd_p_qd(:, :, i);
        p_qd_mid_p_qnext = -(dt/8)*p_qdd_p_q(:, :, i+1);
        p_qd_mid_p_qdnext = eye(nstates/2)/2 - (dt/8)*p_qdd_p_qd(:, :, i+1);
        
        p_qd_mid_p_u = (dt/8)*p_qdd_p_Uall(:, torso_dof+1:end, i);
        p_qd_mid_p_unext = -(dt/8)*p_qdd_p_Uall(:, torso_dof+1:end, i+1);
        
        dyn_gradceq(nstates*(i-1)+1 : nstates*(i-1)+nstates/2, nstates*(i-1)+1 : nstates*(i+1)) = ...
            [-eye(nstates/2) - (4*dt/6)*p_qd_mid_p_q, ...
             -(dt/6)*eye(nstates/2) - (4*dt/6)*p_qd_mid_p_qd, ...
             eye(nstates/2) - (4*dt/6)*p_qd_mid_p_qnext, ...
             -(dt/6)*eye(nstates/2) - (4*dt/6)*p_qd_mid_p_qdnext];

        dyn_gradceq(nstates*(i-1)+1 : nstates*(i-1)+nstates/2, ...
            [nstates*(nsteps+1)+(2*i-2)*ninputs+1 : nstates*(nsteps+1)+(2*i-1)*ninputs, ...
                nstates*(nsteps+1)+(2*i)*ninputs+1 : nstates*(nsteps+1)+(2*i+1)*ninputs]) = [-(4*dt/6)*p_qd_mid_p_u, -(4*dt/6)*p_qd_mid_p_unext];
        
        dyn_gradceq(nstates*(i-1)+nstates/2+1 : nstates*i, nstates*(i-1)+1 : nstates*(i+1)) = ...
            [-(dt/6)*p_qdd_p_q(:,:,i) - (4*dt/6)*(p_qdd_mid_p_q_mid(:,:,i)*p_q_mid_p_q + p_qdd_mid_p_qd_mid(:,:,i)*p_qd_mid_p_q), ...
             -eye(nstates/2) - (dt/6)*p_qdd_p_qd(:,:,i) - (4*dt/6)*(p_qdd_mid_p_q_mid(:,:,i)*p_q_mid_p_qd + p_qdd_mid_p_qd_mid(:,:,i)*p_qd_mid_p_qd), ...
             -(dt/6)*p_qdd_p_q(:,:,i+1) - (4*dt/6)*(p_qdd_mid_p_q_mid(:,:,i)*p_q_mid_p_qnext + p_qdd_mid_p_qd_mid(:,:,i)*p_qd_mid_p_qnext), ...
             eye(nstates/2) - (dt/6)*p_qdd_p_qd(:,:,i+1) - (4*dt/6)*(p_qdd_mid_p_q_mid(:,:,i)*p_q_mid_p_qdnext + p_qdd_mid_p_qd_mid(:,:,i)*p_qd_mid_p_qdnext)];
        
        dyn_gradceq(nstates*(i-1)+nstates/2+1 : nstates*i, nstates*(nsteps+1)+(2*i-2)*ninputs+1 : nstates*(nsteps+1)+(2*i+1)*ninputs) = ...
            [-(dt/6)*p_qdd_p_Uall(:, torso_dof+1:end, i) - (4*dt/6)*p_qdd_mid_p_qd_mid(:,:,i)*p_qd_mid_p_u, ...
             -(4*dt/6)*p_qdd_mid_p_Uall_mid(:, torso_dof+1:end, i), ...
             -(dt/6)*p_qdd_p_Uall(:, torso_dof+1:end, i+1) - (4*dt/6)*p_qdd_mid_p_qd_mid(:,:,i)*p_qd_mid_p_unext];

        switch urdf_filename_drd
            case 'three-link_2JDOF_tail_robot_drd_3DOF_float.urdf'
                coll_c(3*(i-1)+1 : 3*i, 1) = three_link_2JDOF_selfcollision_con_3DOF_float(q(:, i));
                coll_gradc(3*(i-1)+1 : 3*i, nstates*(i-1)+1 : nstates*(i-1)+nstates/2) = three_link_2JDOF_selfcollision_con_grad_3DOF_float(q(:, i));
    
            case 'four-link_2JDOF_tail_robot_drd_3DOF_float.urdf'
                coll_c(6*(i-1)+1 : 6*i, 1) = four_link_2JDOF_selfcollision_con_3DOF_float(q(:, i));
                coll_gradc(6*(i-1)+1 : 6*i, nstates*(i-1)+1 : nstates*(i-1)+nstates/2) = four_link_2JDOF_selfcollision_con_grad_3DOF_float(q(:, i));
            
            case 'five-link_2JDOF_tail_robot_drd_3DOF_float.urdf'
                coll_c(9*(i-1)+1 : 9*i, 1) = five_link_2JDOF_selfcollision_con_3DOF_float(q(:, i));
                coll_gradc(9*(i-1)+1 : 9*i, nstates*(i-1)+1 : nstates*(i-1)+nstates/2) = five_link_2JDOF_selfcollision_con_grad_3DOF_float(q(:, i));
        
            case 'six-link_2JDOF_tail_robot_drd_3DOF_float.urdf'
                coll_c(12*(i-1)+1 : 12*i, 1) = six_link_2JDOF_selfcollision_con_3DOF_float(q(:, i));
                coll_gradc(12*(i-1)+1 : 12*i, nstates*(i-1)+1 : nstates*(i-1)+nstates/2) = six_link_2JDOF_selfcollision_con_grad_3DOF_float(q(:, i));
                
        end
        
        xmid_intergradc1 = [ind_mtx*p_q_mid_p_q, ind_mtx*p_q_mid_p_qd, ind_mtx*p_q_mid_p_qnext, ind_mtx*p_q_mid_p_qdnext;
                           ind_mtx*p_qd_mid_p_q, ind_mtx*p_qd_mid_p_qd, ind_mtx*p_qd_mid_p_qnext, ind_mtx*p_qd_mid_p_qdnext];

        xmid_intergradc2 = [ind_mtx*p_qd_mid_p_u, ind_mtx*p_qd_mid_p_unext];
        
        xmid_gradc(ninputs*4*(i-1)+1 : ninputs*4*i, nstates*(i-1)+1 : nstates*(i+1)) = [xmid_intergradc1; -xmid_intergradc1];

        xmid_gradc([ninputs*(4*i-3)+1 : ninputs*(4*i-2),  ninputs*(4*i-1)+1 : ninputs*4*i], ...
            [nstates*(nsteps+1)+(2*i-2)*ninputs+1 : nstates*(nsteps+1)+(2*i-1)*ninputs, ...
                    nstates*(nsteps+1)+(2*i)*ninputs+1 : nstates*(nsteps+1)+(2*i+1)*ninputs]) = [xmid_intergradc2; -xmid_intergradc2];

    end

    pow_c = sum((Uact').^2, 2) - power_limit*ones(2*nsteps+1, 1);

    Uact_cell = num2cell(Uact', 2);
    pow_gradc(:, (nsteps+1)*nstates+1 : end) = 2*blkdiag(Uact_cell{:});

    %%
    c = [coeff.pow * pow_c;
         coeff.coll * coll_c;
         coeff.xmid * xmid_c];

    gradc = [coeff.pow * pow_gradc;
             coeff.coll * coll_gradc;
             coeff.xmid * xmid_gradc]';

    ceq = coeff.dyn * dyn_ceq;
    gradceq = coeff.dyn * dyn_gradceq';

end

    