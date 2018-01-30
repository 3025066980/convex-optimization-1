%problem 2 of Linear Algebra & Convex term project
randn('state',1);
m=3;
n=2;
alpha=0.01;
beta=0.5;
MAXITERS=1000;
stop=1e-3;
% generate random problem
A=randn(m,n);
b=randn(m,1);
% gradient method
vals=[];
x=zeros(n,1);
steps=[x];
for iter=1:MAXITERS
    val=x'*A'*A*x-2*x'*A'*b+b'*b;
    vals=[vals val];
    grad=2*A'*A*x-2*A'*b;
    v=-grad;
    if norm(grad,2)<stop
        break; 
    end;
    t=1;
    temp=(x+t*v)'*A'*A*(x+t*v)-2*(x+t*v)'*A'*b+b'*b;
    while(temp>(val+alpha*t*grad'*grad))
        t=beta*t;
        temp=(x+t*v)'*A'*A*(x+t*v)-2*(x+t*v)'*A'*b+b'*b;
    end;
    x=x+t*v;
    steps=[steps,x];
end;
figure;
[row col]=size(vals);
optimal=vals(1,col);
plot(1:col,vals(1,1:col)-optimal);
title('f(x)-optimal vs iteration');
figure;
plot(steps(1,:),steps(2,:));
hold;
%cvx prove
cvx_begin quiet
  variable answer(2)
  minimize(answer'*A'*A*answer-2*answer'*A'*b+b'*b)
cvx_end
%contour plot
syms a;syms b
B=A'*A;C=2*A'*b;c=b'*b;
f=a^2*B(1,1)+a*b*(B(1,2)+B(2,1))+b^2*B(2,2)-a*C(1)-b*C(2)+c;
domain=[min(steps(1,:))-0.2,max(steps(1,:))+0.2,min(steps(2,:))-0.2,max(steps(2,:))+0.2];
ezcontour(f,domain);
title('Iterations');
hold;