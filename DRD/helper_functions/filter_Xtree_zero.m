function model = filter_Xtree_zero(model, threshold)

% Due to rounding error, there could be some number in Xtree that are very
% close to 0, -1, 1.
% We directly round them to 0, -1, 1 so that the symbolic computation could
% be easier.

if nargin == 1
    threshold = 1e-7;
end

if isa(model, 'double')
    R = model(1:3,1:3);
    p = model(4:6,1:3);

    R(abs(R) <= threshold) = 0;
    R(abs(R - 1) <= threshold) = 1;
    R(abs(R + 1) <= threshold) = -1;
    p(abs(p) <= threshold) = 0;

    model(1:3,1:3) = R;
    model(4:6,4:6) = R;
    model(4:6,1:3) = p;
else
    for i = 1:model.NB
        Xtree = model.Xtree{i};
        R = Xtree(1:3,1:3);
        p = Xtree(4:6,1:3);
    
        R(abs(R) <= threshold) = 0;
        R(abs(R - 1) <= threshold) = 1;
        R(abs(R + 1) <= threshold) = -1;
        p(abs(p) <= threshold) = 0;
    
        model.Xtree{i}(1:3,1:3) = R;
        model.Xtree{i}(4:6,4:6) = R;
        model.Xtree{i}(4:6,1:3) = p;
    end
end