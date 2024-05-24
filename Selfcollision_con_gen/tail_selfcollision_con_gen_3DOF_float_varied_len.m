clear; close all;

%%
urdf_subdir = '3DOF_floatingbase';
urdf_sets = {'three-link_2JDOF_tail_robot_robtool_3DOF_float.urdf',...
             'four-link_2JDOF_tail_robot_robtool_3DOF_float.urdf', ...
             'five-link_2JDOF_tail_robot_robtool_3DOF_float.urdf', ...
             'six-link_2JDOF_tail_robot_robtool_3DOF_float.urdf'};

% sphere parameters
torso_sph_r = sqrt(2)/6;
ee_sph_r = sqrt(2)/20;

torso_len = 1;

for i = 1 : length(urdf_sets)
    urdf_filename = urdf_sets{i};            

    %%
    robot = importrobot(['..', filesep, 'URDF', filesep, 'urdf_for_robtool', filesep, urdf_subdir, filesep, urdf_filename]);
    robot.DataFormat = 'column';
    
    %%
    torso_dof = 3;
    
    nstates = 2 * length(robot.homeConfiguration);
    ninputs = length(robot.homeConfiguration) - torso_dof;
    
    %%
    q = sym('q', [nstates/2, 1], 'real');
    
    torso_sph1_ctr = eul2rotm_func(q(1:3))*([0; -torso_len/3; 0]);
    torso_sph2_ctr = eul2rotm_func(q(1:3))*([0; 0; 0]);
    torso_sph3_ctr = eul2rotm_func(q(1:3))*([0; torso_len/3; 0]);
    
    %%
    switch urdf_filename
        case {'three-link_2JDOF_tail_robot_robtool_3DOF_float.urdf'}
            link_num = 3;
            jdof = 2;
        case {'four-link_2JDOF_tail_robot_robtool_3DOF_float.urdf'}
            link_num = 4;
            jdof = 2;
        case {'five-link_2JDOF_tail_robot_robtool_3DOF_float.urdf'}
            link_num = 5;
            jdof = 2;
        case {'six-link_2JDOF_tail_robot_robtool_3DOF_float.urdf'}
            link_num = 6;
            jdof = 2;
    end

    %%
    link_len = sym('l', [link_num, 1], 'real');
    
    %%
    if link_num == 3
        joint1_pos = eul2rotm_func(q(1:3))*([0; torso_len/2; 0]);
        joint2_pos = joint1_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*([0; link_len(1); 0]);
        joint3_pos = joint2_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*([0; link_len(2); 0]);

        ee_sph_ctr = joint3_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*...
                     eul2rotm_func([q(8); sym(0); q(9)])*([0; link_len(3); 0]);

    elseif link_num == 4
        joint1_pos = eul2rotm_func(q(1:3))*([0; torso_len/2; 0]);
        joint2_pos = joint1_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*([0; link_len(1); 0]);
        joint3_pos = joint2_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*([0; link_len(2); 0]);
        joint4_pos = joint3_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*...
                     eul2rotm_func([q(8); sym(0); q(9)])*([0; link_len(3); 0]);
        
        ee_sph_ctr = joint4_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*...
                     eul2rotm_func([q(8); sym(0); q(9)])*eul2rotm_func([q(10); sym(0); q(11)])*([0; link_len(4); 0]);

    elseif link_num == 5
        joint1_pos = eul2rotm_func(q(1:3))*([0; torso_len/2; 0]);
        joint2_pos = joint1_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*([0; link_len(1); 0]);
        joint3_pos = joint2_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*([0; link_len(2); 0]);
        joint4_pos = joint3_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*...
                     eul2rotm_func([q(8); sym(0); q(9)])*([0; link_len(3); 0]);
        joint5_pos = joint4_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*...
                     eul2rotm_func([q(8); sym(0); q(9)])*eul2rotm_func([q(10); sym(0); q(11)])*([0; link_len(4); 0]);
        
        ee_sph_ctr = joint5_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*...
                     eul2rotm_func([q(8); sym(0); q(9)])*eul2rotm_func([q(10); sym(0); q(11)])*eul2rotm_func([q(12); sym(0); q(13)])*([0; link_len(5); 0]);
    
    elseif link_num == 6
        joint1_pos = eul2rotm_func(q(1:3))*([0; torso_len/2; 0]);
        joint2_pos = joint1_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*([0; link_len(1); 0]);
        joint3_pos = joint2_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*([0; link_len(2); 0]);
        joint4_pos = joint3_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*...
                     eul2rotm_func([q(8); sym(0); q(9)])*([0; link_len(3); 0]);
        joint5_pos = joint4_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*...
                     eul2rotm_func([q(8); sym(0); q(9)])*eul2rotm_func([q(10); sym(0); q(11)])*([0; link_len(4); 0]);
        joint6_pos = joint5_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*...
                     eul2rotm_func([q(8); sym(0); q(9)])*eul2rotm_func([q(10); sym(0); q(11)])*eul2rotm_func([q(12); sym(0); q(13)])*([0; link_len(5); 0]);

        ee_sph_ctr = joint6_pos + eul2rotm_func(q(1:3))*eul2rotm_func([q(4); sym(0); q(5)])*eul2rotm_func([q(6); sym(0); q(7)])*...
                     eul2rotm_func([q(8); sym(0); q(9)])*eul2rotm_func([q(10); sym(0); q(11)])*eul2rotm_func([q(12); sym(0); q(13)])*...
                     eul2rotm_func([q(14); sym(0); q(15)])*([0; link_len(6); 0]);
    
    end
    
    %%

    if link_num == 3
        selfcoll_con(1) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph1_ctr)'*(ee_sph_ctr - torso_sph1_ctr);
        selfcoll_con(2) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph2_ctr)'*(ee_sph_ctr - torso_sph2_ctr);
        selfcoll_con(3) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph3_ctr)'*(ee_sph_ctr - torso_sph3_ctr);
    
    elseif link_num == 4
        selfcoll_con(1) = (torso_sph_r + ee_sph_r)^2 - (joint4_pos - torso_sph1_ctr)'*(joint4_pos - torso_sph1_ctr);
        selfcoll_con(2) = (torso_sph_r + ee_sph_r)^2 - (joint4_pos - torso_sph2_ctr)'*(joint4_pos - torso_sph2_ctr);
        selfcoll_con(3) = (torso_sph_r + ee_sph_r)^2 - (joint4_pos - torso_sph3_ctr)'*(joint4_pos - torso_sph3_ctr);

        selfcoll_con(4) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph1_ctr)'*(ee_sph_ctr - torso_sph1_ctr);
        selfcoll_con(5) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph2_ctr)'*(ee_sph_ctr - torso_sph2_ctr);
        selfcoll_con(6) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph3_ctr)'*(ee_sph_ctr - torso_sph3_ctr);
    
    elseif link_num == 5
        selfcoll_con(1) = (torso_sph_r + ee_sph_r)^2 - (joint4_pos - torso_sph1_ctr)'*(joint4_pos - torso_sph1_ctr);
        selfcoll_con(2) = (torso_sph_r + ee_sph_r)^2 - (joint4_pos - torso_sph2_ctr)'*(joint4_pos - torso_sph2_ctr);
        selfcoll_con(3) = (torso_sph_r + ee_sph_r)^2 - (joint4_pos - torso_sph3_ctr)'*(joint4_pos - torso_sph3_ctr);
    
        selfcoll_con(4) = (torso_sph_r + ee_sph_r)^2 - (joint5_pos - torso_sph1_ctr)'*(joint5_pos - torso_sph1_ctr);
        selfcoll_con(5) = (torso_sph_r + ee_sph_r)^2 - (joint5_pos - torso_sph2_ctr)'*(joint5_pos - torso_sph2_ctr);
        selfcoll_con(6) = (torso_sph_r + ee_sph_r)^2 - (joint5_pos - torso_sph3_ctr)'*(joint5_pos - torso_sph3_ctr);

        selfcoll_con(7) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph1_ctr)'*(ee_sph_ctr - torso_sph1_ctr);
        selfcoll_con(8) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph2_ctr)'*(ee_sph_ctr - torso_sph2_ctr);
        selfcoll_con(9) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph3_ctr)'*(ee_sph_ctr - torso_sph3_ctr);

    elseif link_num == 6
        selfcoll_con(1) = (torso_sph_r + ee_sph_r)^2 - (joint4_pos - torso_sph1_ctr)'*(joint4_pos - torso_sph1_ctr);
        selfcoll_con(2) = (torso_sph_r + ee_sph_r)^2 - (joint4_pos - torso_sph2_ctr)'*(joint4_pos - torso_sph2_ctr);
        selfcoll_con(3) = (torso_sph_r + ee_sph_r)^2 - (joint4_pos - torso_sph3_ctr)'*(joint4_pos - torso_sph3_ctr);
    
        selfcoll_con(4) = (torso_sph_r + ee_sph_r)^2 - (joint5_pos - torso_sph1_ctr)'*(joint5_pos - torso_sph1_ctr);
        selfcoll_con(5) = (torso_sph_r + ee_sph_r)^2 - (joint5_pos - torso_sph2_ctr)'*(joint5_pos - torso_sph2_ctr);
        selfcoll_con(6) = (torso_sph_r + ee_sph_r)^2 - (joint5_pos - torso_sph3_ctr)'*(joint5_pos - torso_sph3_ctr);

        selfcoll_con(7) = (torso_sph_r + ee_sph_r)^2 - (joint6_pos - torso_sph1_ctr)'*(joint6_pos - torso_sph1_ctr);
        selfcoll_con(8) = (torso_sph_r + ee_sph_r)^2 - (joint6_pos - torso_sph2_ctr)'*(joint6_pos - torso_sph2_ctr);
        selfcoll_con(9) = (torso_sph_r + ee_sph_r)^2 - (joint6_pos - torso_sph3_ctr)'*(joint6_pos - torso_sph3_ctr);

        selfcoll_con(10) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph1_ctr)'*(ee_sph_ctr - torso_sph1_ctr);
        selfcoll_con(11) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph2_ctr)'*(ee_sph_ctr - torso_sph2_ctr);
        selfcoll_con(12) = (torso_sph_r + ee_sph_r)^2 - (ee_sph_ctr - torso_sph3_ctr)'*(ee_sph_ctr - torso_sph3_ctr);
    
    end
    
    pcon_pq = jacobian(selfcoll_con.', q);
    pcon_plen = jacobian(selfcoll_con.', link_len);

    %%
    if link_num == 3
        matlabFunction(selfcoll_con.', 'Vars', {q, link_len}, 'File', ['selfcollision_con_varied_len', filesep, 'three_link_2JDOF_selfcollision_con_3DOF_float_varied_len', ]);
        matlabFunction(pcon_pq, pcon_plen, 'Vars', {q, link_len}, 'File', ['selfcollision_con_varied_len', filesep, 'three_link_2JDOF_selfcollision_con_grad_3DOF_float_varied_len']);
    
    elseif link_num == 4
        matlabFunction(selfcoll_con.', 'Vars', {q, link_len}, 'File', ['selfcollision_con_varied_len', filesep, 'four_link_2JDOF_selfcollision_con_3DOF_float_varied_len', ]);
        matlabFunction(pcon_pq, pcon_plen, 'Vars', {q, link_len}, 'File', ['selfcollision_con_varied_len', filesep, 'four_link_2JDOF_selfcollision_con_grad_3DOF_float_varied_len']);
    
    elseif link_num == 5
        matlabFunction(selfcoll_con.', 'Vars', {q, link_len}, 'File', ['selfcollision_con_varied_len', filesep, 'five_link_2JDOF_selfcollision_con_3DOF_float_varied_len', ]);
        matlabFunction(pcon_pq, pcon_plen, 'Vars', {q, link_len}, 'File', ['selfcollision_con_varied_len', filesep, 'five_link_2JDOF_selfcollision_con_grad_3DOF_float_varied_len']);
    
    elseif link_num == 6
        matlabFunction(selfcoll_con.', 'Vars', {q, link_len}, 'File', ['selfcollision_con_varied_len', filesep, 'six_link_2JDOF_selfcollision_con_3DOF_float_varied_len', ]);
        matlabFunction(pcon_pq, pcon_plen, 'Vars', {q, link_len}, 'File', ['selfcollision_con_varied_len', filesep, 'six_link_2JDOF_selfcollision_con_grad_3DOF_float_varied_len']);

    end
       
end

%%

