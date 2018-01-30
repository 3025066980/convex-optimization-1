%problem 1 of Linear Algebra & Convex term project
%n=50;
%n=100;
n=500;
S=5;        %number of nonzero elements in x
initial=3;  %initial value for k
increment=50; %increment for k
rounds=20;    %rounds for Monte Carlo approach
start=43;     %start position of nonzero elements in x
P=zeros(increment,1);   %probability vector
c=1;
for k=initial:(initial+increment)
    j=0;
    for i=1:rounds
        x=zeros(n,1);
        x(start:start-1+S,1)=rand(S,1);       %give values to S nonzero elements
        A=2*binornd(1,0.5,k,n)-ones(k,n);     %generate matrix A
        y=A*x;                                %compute y in the constraint
        %recovery of x
        cvx_begin quiet
           variable d(n)
           minimize(norm(d,1))
           subject to
               y==A*d    
        cvx_end
        %compare the recovered x to the original x
        if abs(sum(d-x))<0.01
            j=j+1;
        else
            j=j;
        end
    end
    P(c)=j/rounds;    %compute the perfect recovery probability
    c=c+1;
end
plot(initial:(initial+increment),P);
title('Probability of perfect recovery versus k,n=500');