function [robot] = removeFixedJoints(robot)
% Update inertias and transforms
for i = 1:robot.NB
    if ~strcmp(robot.jtype{i}, 'fixed')
        continue;
    end
    
    parent = robot.parent(i);
    children = find(robot.parent == i);
    
    % Update the transforms of the children
    for j = children
        robot.parent(j) = parent;
        robot.Xtree{j} = robot.Xtree{j}*robot.Xtree{i};
    end
    
    % Absorb the body into the parent body if it exists
    if parent > 0
        % Add inertia to parent
        X = robot.Xtree{i};
        I = X'*robot.I{i}*X;
        robot.I{parent} = robot.I{parent} + I;
    end
end

% Remove the joints
i = 1;
while i <= robot.NB
    if ~strcmp(robot.jtype{i}, 'fixed')
        i = i + 1;
    else
        robot.NB = robot.NB - 1;
        robot.parent(i) = [];
        robot.jointNames(i) = [];
        robot.Xtree(i) = [];
        robot.I(i) = [];
        robot.jtype(i) = [];
        
        if isfield(robot, 'transmissionInertia')
            robot.transmissionInertia(i) = [];
        end
        if isfield(robot, 'damping')
            robot.damping(i) = [];
        end
        if isfield(robot, 'friction')
            robot.friction(i) = [];
        end
        if isfield(robot, 'upperLim')
            robot.upperLim(i) = [];
            robot.lowerLim(i) = [];
        end
        
        for j = i:robot.NB
            if robot.parent(j) >= i
                robot.parent(j) = robot.parent(j) - 1;
            end
        end
    end
end
end