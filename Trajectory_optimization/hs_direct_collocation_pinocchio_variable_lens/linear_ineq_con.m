function [A, b] = linear_ineq_con(nsteps, nstates, ninputs, link_num, dt, inputrate_range, inputacc_range)

    ndec = (nsteps + 1)*nstates + (nsteps*2 + 1)*ninputs + link_num;

    inputrate_lb = inputrate_range(:, 1);
    inputrate_ub = inputrate_range(:, 2);

    if ~isempty(inputacc_range)
        inputacc_lb = inputacc_range(:, 1);
        inputacc_ub = inputacc_range(:, 2);
    end

    %%
    A_inputrate = sparse((nsteps*2)*ninputs*2, ndec);
    b_inputrate = sparse((nsteps*2)*ninputs*2, 1);

    if ~isempty(inputacc_range)
        A_inputacc = sparse((nsteps*2-1)*ninputs*2, ndec);
        b_inputacc = sparse((nsteps*2-1)*ninputs*2, 1);
    end

    for i = 1 : nsteps*2
        A_inputrate(ninputs*(2*i-2)+1 : ninputs*(2*i-1), (nsteps+1)*nstates+ninputs*(i-1)+1 : (nsteps+1)*nstates+ninputs*(i+1)) = ...
            [-eye(ninputs), eye(ninputs)];
        b_inputrate(ninputs*(2*i-2)+1 : ninputs*(2*i-1), 1) = inputrate_ub.*dt/2;
        
        A_inputrate(ninputs*(2*i-1)+1 : ninputs*(2*i), (nsteps+1)*nstates+ninputs*(i-1)+1 : (nsteps+1)*nstates+ninputs*(i+1)) = ...
            [eye(ninputs), -eye(ninputs)];
        b_inputrate(ninputs*(2*i-1)+1 : ninputs*(2*i), 1) = -inputrate_lb.*dt/2;

        if ~isempty(inputacc_range)

            if i <= nsteps*2 - 1
                A_inputacc(ninputs*(2*i-2)+1 : ninputs*(2*i-1), (nsteps+1)*nstates+ninputs*(i-1)+1 : (nsteps+1)*nstates+ninputs*(i+2)) = ...
                    [eye(ninputs), -2*eye(ninputs), eye(ninputs)];
                b_inputacc(ninputs*(2*i-2)+1 : ninputs*(2*i-1), 1) = inputacc_ub.*(dt/2)^2;
    
                A_inputacc(ninputs*(2*i-1)+1 : ninputs*(2*i), (nsteps+1)*nstates+ninputs*(i-1)+1 : (nsteps+1)*nstates+ninputs*(i+2)) = ...
                    [-eye(ninputs), 2*eye(ninputs), -eye(ninputs)];
                b_inputacc(ninputs*(2*i-1)+1 : ninputs*(2*i), 1) = -inputacc_lb.*(dt/2)^2;
    
            end
        end              
    
    end
        
    %%
    A = [A_inputrate; A_inputacc]; 
    b = [b_inputrate; b_inputacc];

end

