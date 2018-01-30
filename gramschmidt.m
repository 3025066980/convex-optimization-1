%linear algebra & convex opt homework1
%Find an orthonormal for the range of a matrix,and the null space of its
%transpose
function [RAbasis,RATbasis]=gramschmidt(a)
[row,column]=size(a);
column2=rank(a);
RAbasis=zeros(row,column2);
for i=1:column2
    temp=0;
    for j=1:i
        temp=RAbasis(:,j)'*a(:,i)*RAbasis(:,j);
    end
    RAbasis(:,i)=a(:,i)-temp;
    RAbasis(:,i)=RAbasis(:,i)/norm(RAbasis(:,i));
end
if column2>row
    RATbasis=0;
else
    RATbasis=linsolve(RAbasis',zeros(column2,1));
end
orth(a)    %basis of range of a
null(a')   %basis of nullspace of a transpose
