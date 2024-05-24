clear; clc;

%% set up the search path
addpath(genpath(['..', filesep, 'DRD']));

%% specify the urdf
% the direcotry in which the urdf file is stored
urdf_subdir = '3DOF_floatingbase';
% the name of the urdf file
urdf_filename = 'two-link_2JDOF_tail_robot_robtool_3DOF_float.urdf';
% the number of tail links (should match the urdf name)
num_links = 2;

%% import urdf to create a model
model = create_model_from_urdf(['..', filesep, 'URDF', filesep, 'urdf_for_robtool', filesep, urdf_subdir, filesep, urdf_filename]);
model.gravity = [0; 0; 0];

if isfield(model, 'transmissionInertia')
    model = rmfield(model, 'transmissionInertia');
end
if isfield(model, 'damping')
    model = rmfield(model, 'damping');
end
if isfield(model, 'friction')
    model = rmfield(model, 'friction');
end

%% generate symbolic variables for system states
% the total length of the tail is 1.5; the total mass of the the tail is 3
% rho represents the ratio: mass/length
rho = 2;

model.NB = 3 + 2 * num_links;

link_lengths = sym('l', [num_links, 1], 'real');

q = sym('q', [model.NB, 1], 'real');
qd = sym('qd', [model.NB, 1], 'real');
qdd = sym('qdd', [model.NB, 1], 'real');

%% symbolize the model parameters
for i = 1 : num_links
    % tail link i virtual
    model.parent(3 + 2*i - 1) = 3 + 2*i - 2;
    model.jtype{3 + 2*i - 1} = 'Rx';
    model.I{3 + 2*i - 1} = zeros(6);
    if i > 1
        model.Xtree{3 + 2*i - 1} = plux(eye(3), [0; link_lengths(i); 0]);
    end

    % tail link i
    model.parent(3 + 2*i) = 3 + 2*i - 1;
    model.jtype{3 + 2*i} = 'Rz';
    
    mass = rho * link_lengths(i);
    
    com = [0, 0.5*link_lengths(i), 0];
    box_dim = [0.1, link_lengths(i), 0.1];

    inertia = [mass*(box_dim(2)^2 + box_dim(3)^2)/12, 0, 0; ...
               0, mass*(box_dim(1)^2 + box_dim(3)^2)/12, 0; ...
               0, 0, mass*(box_dim(1)^2 + box_dim(2)^2)/12];
    
    model.I{3 + 2*i} = mcI(mass, com, inertia);
    model.Xtree{3 + 2*i} = eye(6);

end

%% set up the prefix of the function to be generated
switch urdf_filename
    case 'one-link_2JDOF_tail_robot_robtool_3DOF_float.urdf'
        funcname2save_pref = 'one_link_2JDOF_varied_len';
        
    case 'two-link_2JDOF_tail_robot_robtool_3DOF_float.urdf'
        funcname2save_pref = 'two_link_2JDOF_varied_len';
        
    case 'three-link_2JDOF_tail_robot_robtool_3DOF_float.urdf'
        funcname2save_pref = 'three_link_2JDOF_varied_len';

    case 'four-link_2JDOF_tail_robot_robtool_3DOF_float.urdf'
        funcname2save_pref = 'four_link_2JDOF_varied_len';
    
    case 'five-link_2JDOF_tail_robot_robtool_3DOF_float.urdf'
        funcname2save_pref = 'five_link_2JDOF_varied_len';
    
    otherwise
        error('Undefined urdf name.')

end

%% set up the directory that stores the generated functions
ptau_plen_savedir = 'ptau_plen';

if ~isfolder(ptau_plen_savedir)
    mkdir(ptau_plen_savedir);
end

%% generate and export the jacobian
% calculate the symbolic system control inputs using inverse dynamics
tau = ID(model, q, qd, qdd);

% calculate the jacobian
ptau_plen = jacobian(tau, link_lengths);

% export the jacobian
matlabFunction(ptau_plen, 'File', [ptau_plen_savedir, filesep, funcname2save_pref, '_ptau_plen.m'], 'Vars', {q, qd, qdd, link_lengths});

%%

