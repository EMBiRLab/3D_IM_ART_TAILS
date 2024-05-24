function [lb, ub] = bound_con(nsteps, state_range, input_range)

    lbNub = [repmat(state_range, nsteps+1, 1); repmat(input_range, nsteps*2 + 1, 1)];
    lb = lbNub(:, 1);
    ub = lbNub(:, 2);

end
