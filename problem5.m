%Linear Algebra & Convex term project problem 5
%A=rand(10,2);
%B=rand(10,2);
A=[0.2 0.5;0.3 0.7;1 2;1 3;0.9 1.4;0.4 0.3;1.1 1.2;0.8 2.2];
B=[-1.1 2;-0.3 0.2;-0.7 0.9;-2 0.4;-1.3 1.4;-2.2 2.4;-0.7 1.7;-1.4 .6];
K=convhull(B);             %generate the smallest convex polyhedra that covers all the points in set B
len=length(K)-1;
hull=zeros(len,2);
%compute the slopes and intercepts of the boundary and form matrix hull
for i=1:len-1
    slope=(B(K(i),2)-B(K(i+1),2))/(B(K(i),1)-B(K(i+1),1));
    intercept=B(K(i),2)-slope*B(K(i),1);
    hull(i,1:2)=[slope intercept];
end
slope=(B(K(len),2)-B(K(1),2))/(B(K(len),1)-B(K(1),1));
intercept=B(K(len),2)-slope*B(K(len),2);
hull(len,1:2)=[slope intercept];
%find a point in B that is not on the boundary
R=sortrows(K);
j=2;
while 1
    if (R(j)+1)~=R(j+1)
        break
    else
        j=j+1;
    end
end
point=B(R(j)+1,1:2);
%compute the values to determine whether a point is inside the polyhedra
compare1=ones(len,1)*point(2);
compare2=point(1)*hull(:,1)+hull(:,2);
compare=compare2-compare1;
sum=0;
%check every point in set A
for i=1:length(A)
    compare3=A(i,1)*hull(:,1)+hull(:,2)-ones(len,1)*A(i,2);
    acc=0;
    for j=1:len
        if compare(j)>0 & compare3(j)>0
            acc=acc+1;
        elseif compare(j)<0 & compare3(j)<0
            acc=acc+1;
        else
            acc=acc;
        end
    end
    if acc==len
        sum=sum+1;
    end
end
if sum>0
    disp('Infeasible');
    plot(A(:,1),A(:,2),'ro');
    hold;
    plot(B(:,1),B(:,2),'ko');
    plot(B(K,1),B(K,2),'k');
    title('Infeasible');
    hold;
else
    disp('Feasible');
    plot(A(:,1),A(:,2),'ro');
    hold;
    plot(B(:,1),B(:,2),'ko');
    plot(B(K,1),B(K,2),'k');
    title('Feasible');
    hold;
end
