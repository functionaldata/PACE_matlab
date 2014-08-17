%  Example for nonlindyn.m

%  Fixed, regular and dense design, 50 subjects with 51 measurements with equal distance  from each subject.
%  Time interval is [0,10].
%  The true dynamic is y'(t)^2=1-y(t)^2
n=50;
t=cell(1,n);
y=cell(1,n);
ti=0:0.2:10;

for i=1:n
  t{i}=ti;
  y{i}=sin(t{i})+0.1*randn(1,1);
end

[f,dyn_grid]=nonlindyn(t,y,'gauss');  