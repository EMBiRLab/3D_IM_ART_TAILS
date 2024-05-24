function robot = create_model_from_urdf(urdf_file_name, base_rot_axis, if_floating_base, uncertainty_mcI)

if nargin < 3
    if_floating_base = false;
end

if nargin < 2
    base_rot_axis = 'Rx';
    disable_base_rot_axis = true;
else
    if strcmp(base_rot_axis, 'x') || strcmp(base_rot_axis, 'y') || strcmp(base_rot_axis, 'z')
        base_rot_axis = ['R',base_rot_axis];
        disable_base_rot_axis = false;
    else
        disable_base_rot_axis = true;
    end
end

[robot.name, links, joints, transmissions] = ros_load_urdf(urdf_file_name);

robot.NB = length(links);

robot.parent = 0:(robot.NB-1);
%More sophisticated parent determiniation
for i = 2:robot.NB
    name = links(i).Name;
    for j = 1:robot.NB - 1
        if strcmp(joints(j).Child, name)
            break
        end
    end
    parent = joints(j).Parent;
    for j = 1:robot.NB
        if strcmp(links(j).Name, parent)
            break
        end
    end
    robot.parent(i) = j;
end

for i = 1:robot.NB
    robot.jointNames(i) = string(links(i).Name); %Joints are named the same as links in the MuJoCo model
end

robot.gravity = [0;0;-9.81];

robot.jtype = cell(1,robot.NB);
robot.Xtree = cell(1,robot.NB);
robot.transmissionInertia = zeros(1,robot.NB);
robot.friction = zeros(1,robot.NB);
robot.damping = zeros(1,robot.NB);
robot.I = cell(1,robot.NB);
robot.appearance.body = cell(1,robot.NB);

for i = 1:robot.NB
    if i == 1
        robot.jtype{i} = base_rot_axis;
        robot.Xtree{i} = eye(6);
        robot.transmissionInertia(i) = 0;
        robot.friction(i) = 0;
        robot.damping(i) = 0;
    else
        if ~isempty(joints(i-1).Axis)
            if strcmp(joints(i-1).Type, 'revolute') || strcmp(joints(i-1).Type, 'continuous')
                jtype = 'R';
            else
                jtype = 'P';
            end

            if joints(i-1).Axis(1) ~= 0
                if joints(i-1).Axis(1) == -1
                    robot.jtype{i} = ['-',jtype,'x'];
                else
                    robot.jtype{i} = [jtype,'x'];
                end
            else
                if joints(i-1).Axis(2) ~= 0
                    if joints(i-1).Axis(2) == -1
                        robot.jtype{i} = ['-',jtype,'y'];
                    else
                        robot.jtype{i} = [jtype,'y'];
                    end
                else
                    if joints(i-1).Axis(3) ~= 0
                        if joints(i-1).Axis(3) == -1
                            robot.jtype{i} = ['-',jtype,'z'];
                        else
                            robot.jtype{i} = [jtype,'z'];
                        end
                    else
                        robot.jtype{i} = 'fixed';
                    end
                end
            end
        else
            robot.jtype{i} = 'fixed';
        end
        % friction and dampings
        robot.friction(i) = joints(i-1).friction;
        robot.damping(i) = joints(i-1).damping;

        if size(joints(i-1).R,1) == 3 && size(joints(i-1).R,2) == 3
            robot.Xtree{i} = inv(pluho([joints(i-1).R, joints(i-1).Offset'; 0 0 0 1]));
        else
%             robot.Xtree{i} = rotx(joints(i-1).R(1)) * roty(joints(i-1).R(2)) * rotz(joints(i-1).R(3)) * xlt(joints(i-1).Offset);
            robot.Xtree{i} = rotx_clean(joints(i-1).R(1)) * roty_clean(joints(i-1).R(2)) * rotz_clean(joints(i-1).R(3)) * xlt(joints(i-1).Offset);
        end
        % transmission inertia (armature in MJCF)
        robot.transmissionInertia(i) = 0;
        for j = 1:length(transmissions)
            if strcmp(joints(i-1).Name, transmissions(j).Joint)
                robot.transmissionInertia(i) = transmissions(j).Inertia * (transmissions(j).MechanicalReduction)^2;
                break;
            end
        end
    end

    if nargin > 3
        robot.I{i} = mcI_interval(links(i).Mass, links(i).Offset, links(i).Inertia, uncertainty_mcI.um, uncertainty_mcI.uc, uncertainty_mcI.uI);
    else
        robot.I{i} = mcI(links(i).Mass, links(i).Offset, links(i).Inertia);
    end

    robot.appearance.body{i} = { 'cyl', [0 0 0; 0 0.4 0], 0.01 };
end

if if_floating_base
    robot = floatbase(robot);
    robot = removeFixedJoints(robot); % Remove extra fixed joints
    robot = filter_Xtree_zero(robot);
    return;
end

if disable_base_rot_axis
    %Cut out Rx base joint
    robot.parent = robot.parent(2:robot.NB) - 1;
    robot.jointNames = robot.jointNames(2:robot.NB);
    robot.jtype = {robot.jtype{2:robot.NB}};
    robot.Xtree = {robot.Xtree{2:robot.NB}};
    robot.transmissionInertia = robot.transmissionInertia(2:robot.NB);
    robot.friction = robot.friction(2:robot.NB);
    robot.damping = robot.damping(2:robot.NB);
    robot.I = {robot.I{2:robot.NB}};
    robot.appearance.body = {robot.appearance.body{2:robot.NB}};
    robot.NB = robot.NB - 1;
end

robot = removeFixedJoints(robot); % Remove extra fixed joints
robot = filter_Xtree_zero(robot);
