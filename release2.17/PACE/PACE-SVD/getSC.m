function [res invalid regular]=getSC(x, t_x, p_x, y, t_y, p_y, bwccov, nsvd, ngrid, FVE)

if isempty(p_x)
    p_x = setOptions( 'selection_k', nsvd, 'screePlot', 0);
else
    p_x = setVal(p_x, 'selection_k', nsvd);
end
xx = FPCA(x,t_x,p_x);
if isempty(p_y)
    p_y = setOptions('selection_k', nsvd, 'screePlot', 0);
else
    p_y = setVal(p_y, 'selection_k', nsvd);
end

yy = FPCA(y,t_y,p_y);
regular1 = getVal(xx, 'regular');
regular2 = getVal(yy, 'regular');
regular = min(regular1, regular2);
x = getVal(xx,'y');
t_x = getVal(xx,'t');
y = getVal(yy,'y');
t_y = getVal(yy,'t');
if ~isempty(getVal(xx,'no_opt')) && ~isempty(getVal(yy,'no_opt'))
    [res invalid] = getSC1(x, t_x, xx, y, t_y, yy, bwccov, nsvd, regular, ngrid, FVE);
else
    invalid = 1;
    res = [];
end
end