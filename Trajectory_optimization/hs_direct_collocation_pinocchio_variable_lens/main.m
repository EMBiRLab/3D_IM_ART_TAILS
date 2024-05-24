clear; close all;

%%
addpath(genpath(['..', filesep, '..', filesep, 'DRD']));
addpath(genpath(['..', filesep, '..', filesep, 'URDF']));
addpath(genpath(['..', filesep, '..', filesep, 'Selfcollision_con_gen']));
addpath(genpath(['..', filesep, '..', filesep, 'Varied_length_derivative']));
addpath(genpath('mex_files'));

%%
urdf_filename_drd = 'two-link_2JDOF_tail_robot_drd_3DOF_float.urdf';
urdf_filename_pino = strrep(urdf_filename_drd, 'drd', 'pinocchio');

urdf_filepath_drd = which(urdf_filename_drd);
urdf_filepath_pino = which(urdf_filename_pino);

model = create_model_from_urdf(urdf_filepath_drd);
model.gravity = [0; 0; 0];

%% solution dir setup
[~, model_name, ~] = fileparts(urdf_filepath_drd);
sol_dir = ['opt_sol', filesep, model_name];

saveon = true;

%% PRY trajectories
PRY_traj_dir = ['..', filesep, '..', filesep, 'Target_trajectories', filesep, 'T_0_0p002_0p5', filesep];

PRY_traj_name_woidx = 'trial';

traj_num_idx = 1:100;
PRY_traj_set = cell(length(traj_num_idx), 1);

for i = 1 : length(traj_num_idx)
    data_loaded = load([PRY_traj_dir, PRY_traj_name_woidx, num2str(traj_num_idx(i)), '.mat']);
    PRY_traj_set{i, 1} = data_loaded.PRY_traj;
end

%%
t0 = 0; dt = 0.004; tend = 0.5;
T = t0 : dt : tend;

for i = 1 : length(PRY_traj_set)
    t_in_data = PRY_traj_set{i, 1}.t;
    
    if ~isequal(t_in_data, t0 : dt/2 : tend)
        error(['T in the PRY traj', num2str(i) ' does not match the specified T in main.m']);
    end
end

%% state, state velocity, input, and power constraints
torso_dof = 3;
PRY_idx_in_q = [1, 2, 3];

nstates = 2 * model.NB;
ninputs = model.NB - torso_dof;

switch urdf_filename_drd
    case {'one-link_2JDOF_tail_robot_drd_3DOF_float.urdf', 'two-link_2JDOF_tail_robot_drd_3DOF_float.urdf', ...
          'three-link_2JDOF_tail_robot_drd_3DOF_float.urdf', 'four-link_2JDOF_tail_robot_drd_3DOF_float.urdf', ...
          'five-link_2JDOF_tail_robot_drd_3DOF_float.urdf', 'six-link_2JDOF_tail_robot_drd_3DOF_float.urdf'}
        
        njoint = ninputs/2;
        joint_range = [-pi/3, pi/3; -pi/3, pi/3];
        jointspeed_range = [-2*pi, 2*pi; -2*pi, 2*pi];
        state_range = [repmat([-Inf, Inf], torso_dof, 1); repmat(joint_range, njoint, 1); ...
                       repmat([-Inf, Inf], torso_dof, 1); repmat(jointspeed_range, njoint, 1)];
    
    otherwise
        error('Did not find the category to which the urdf belongs.')
end

input_range = repmat([-5, 5], ninputs, 1);
inputrate_range = repmat([-100, 100], ninputs, 1);
inputacc_range = repmat([-500, 500], ninputs, 1);
power_limit = 50;

link_num = njoint;
total_len = 1.5;
linklen_range = repmat([0.2, total_len], link_num, 1);

%%
disp(['Running optimizations for ', model_name, '...']);

opt_coeff.obj = 200;
opt_coeff.con.dyn = 10;
opt_coeff.con.pow = 1e-3;
opt_coeff.con.coll = 1;
opt_coeff.con.xmid = 1;

for i = 1 : length(PRY_traj_set)
    PRY_traj = PRY_traj_set{i, 1};
    PRY_ref = PRY_traj.y;

    PRY_ref_name_widx = [PRY_traj_name_woidx, num2str(traj_num_idx(i))];

    nsteps = length(T) - 1;
    ndec = (nsteps + 1)*nstates + (nsteps*2 + 1)*ninputs + link_num;
    
    x0 = zeros(nstates, 1);
    x0(PRY_idx_in_q) = PRY_ref(:, 1);
    
    objfun = @(z) objective_fun(z, nsteps, nstates, dt, PRY_ref, PRY_idx_in_q, opt_coeff.obj);
    [A, b] = linear_ineq_con(nsteps, nstates, ninputs, link_num, dt, inputrate_range, inputacc_range);
    [Aeq, beq] = linear_eq_con(nsteps, nstates, ninputs, x0, link_num, total_len);
    [lb, ub] = bound_con(nsteps, state_range, input_range, linklen_range);

    nonlcon = @(z) nonlinear_con(z, nsteps, nstates, ninputs, torso_dof, dt, power_limit, state_range, link_num, urdf_filename_drd, urdf_filepath_pino, opt_coeff.con);
    
    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, ...
                                      'SpecifyConstraintGradient', true);
    options.MaxIterations = 4000;
    options.CheckGradients = true;
    options.FiniteDifferenceType = "central";
    options.Display = 'iter';
    options.Algorithm = "interior-point";
    options.ConstraintTolerance = opt_coeff.con.dyn * 1e-6;
    switch options.Algorithm
        case 'interior-point'
            options.MaxFunctionEvaluations = 100 * ndec;
    end
    
    % set up different initial conditions here for the optimization if
    % desired
    z0 = [x0; zeros(nsteps*nstates + (nsteps*2+1)*ninputs, 1); total_len/link_num*ones(link_num, 1)];
    [zsol, fval] = fmincon(objfun, z0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    
    X = reshape(zsol(1 : (nsteps + 1)*nstates, 1), nstates, [])';
    U = reshape(zsol((nsteps + 1)*nstates+1 : (nsteps + 1)*nstates + (nsteps*2 + 1)*ninputs, 1), ninputs, [])';
    L = reshape(zsol(end-link_num+1 : end, 1), link_num, [])';
    
    if saveon
        savedata_dir = [sol_dir, filesep, PRY_ref_name_widx];
        if ~isfolder(savedata_dir)
            mkdir(savedata_dir);
        end
        
        save([savedata_dir, filesep, 'traj_opt_saved_data.mat'], ...
             'urdf_filename_drd', 'model', 't0', 'dt', 'tend', 'T', ...
             'state_range', 'input_range', 'inputrate_range', 'inputacc_range', 'power_limit', 'opt_coeff',...
             'PRY_traj', 'PRY_ref_name_widx', 'options', 'fval', 'zsol', 'X', 'U', 'L');

    end

end

disp('Done!');

%%


