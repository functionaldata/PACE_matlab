function [regular] = isregular(t)
    tt = cell2mat(t);
    f = length(tt)/length(unique(tt))/length(t);
    if f==1
        regular = 2;
    elseif f>.75
        regular = 1;
    else
        regular = 0;
    end
end
