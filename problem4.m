%Linear Algebra term project Problem 4
MAXITERS = 100;        %maximum iterations
ALPHA = 0.01;
BETA = 0.5;
RESTOL = 1e-7;         %stop error
x=[1 2]';              %initial value not in the feasible set
b=[0]';
A=[1 -1];
nu=zeros(1,1);
resdls=[];
vals=[];
for i=1:MAXITERS
    %gradient
    grad=[exp(x(1)+2*x(2))+exp(x(1)-2*x(2))-exp(-x(1)) 2*exp(x(1)+2*x(2))-2*exp(x(1)-2*x(2))]';
    gra22=2*exp(x(1)+2*x(2))-2*exp(x(1)-2*x(2));
    gra24=4*exp(x(1)+2*x(2))+4*exp(x(1)-2*x(2));
    func=exp(x(1)+2*x(2))+exp(x(1)-2*x(2))+exp(-x(1));
    vals=[vals func];
    %Hessian matrix
    grad2=[func gra22;gra22 gra24];
    r=[grad+A'*nu;A*x-b]; 
    resdls=[resdls,norm(r)];
    sol=-[grad2 A';A 0]\r;
    %newton step
    Dx=sol(1:2); 
    %primal and dual step
    Dnu=sol(3);
    if(norm(r)<RESTOL)
         break; 
    end;
    t=1;
    xtemp=x+t*Dx;
    gradtemp=[exp(xtemp(1)+2*xtemp(2))+exp(xtemp(1)-2*xtemp(2))-exp(-xtemp(1)) 2*exp(xtemp(1)+2*xtemp(2))-2*exp(xtemp(1)-2*xtemp(2))]';
    %backtracking
    while(norm([gradtemp+A'*(nu+t*Dnu);A*xtemp-b])>(1-ALPHA*t)*norm(r))
        t=BETA*t;
        xtemp=x+t*Dx;
        gradtemp=[exp(xtemp(1)+2*xtemp(2))+exp(xtemp(1)-2*xtemp(2))-exp(-xtemp(1)) 2*exp(xtemp(1)+2*xtemp(2))-2*exp(xtemp(1)-2*xtemp(2))]';
    end;
    %update x
    x=x+t*Dx;
    nu=nu+t*Dnu;
end;
optimal=vals(1,i);
plot(1:i,vals(1:1:i)-optimal);
title('value - P* versus the number of iterations');
%cvx check
cvx_begin quiet
   variable a(2)
   func=exp(a(1)+2*a(2))+exp(a(1)-2*a(2))+exp(-a(1))
   minimize(func)
   subject to
   a(1)==a(2)
cvx_end
    