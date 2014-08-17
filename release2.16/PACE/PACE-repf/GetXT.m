function [x, t_x, insub] = GetXT(i, n, out1, y, t)
    insub = zeros(1,n);
    for j = 1:n
        id = ismember(t{j}(1,:), out1(i));
        if sum(id) > 0
           insub(j) = 1;
        end
    end
    ni = sum(insub);
    x = cell(1, ni);
    t_x = cell(1, ni);
    k=0;
    for j = 1:n
       if insub(j)==1
           k=k+1;
           id = ismember(t{j}(1,:), out1(i));
           x{k} = y{j}(id);
           t_x{k} = t{j}(2, id);
       end
    end


