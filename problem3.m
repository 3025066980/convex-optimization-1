%Linear Algebra term project problem3
x=[1 1]';     %initial value
stop=1e-6;    %stop error
MAXITE=10;    %maximum iteration rounds
alpha=0.3;   
beta=0.7;
vals=[];      %vals of function output 
for i=1:MAXITE
    func=exp(x(1)+2*x(2))+exp(x(1)-2*x(2))+exp(-x(1));  %function
    gra22=2*exp(x(1)+2*x(2))-2*exp(x(1)-2*x(2));
    gra24=4*exp(x(1)+2*x(2))+4*exp(x(1)-2*x(2));
    gra1=[func-2*exp(-x(1));gra22];                     %first order gradient
    gra2=[func gra22;gra22 gra24];                     %Hessians matrix
    newtonstep=-inv(gra2)*gra1;                        %newton step
    decrement=gra1'*inv(gra2)*gra1;                    %newton decrement
    %check stop criterion
    if decrement<2*stop
        break;
    end
    %backtracking line search
    t=1;                                         
    xtemp=x+t*newtonstep;
    val=exp(xtemp(1)+2*xtemp(2))+exp(xtemp(1)-2*xtemp(2))+exp(-xtemp(1));
    while(val>(func+alpha*t*gra1'*newtonstep))
        t=beta*t;
        xtemp=x+t*newtonstep;
        val=exp(xtemp(1)+2*xtemp(2))+exp(xtemp(1)-2*xtemp(2))+exp(-xtemp(1));
    end
    %update x
    x=x+t*newtonstep;
    vals=[vals func];
end
optimal=vals(1,i-1);                %optimal value
plot(1:i-1,vals(1:1:i-1)-optimal);
title('value - P* versus the number of iterations');
%cvx check
cvx_begin quiet
   variable a(2)
   func=exp(a(1)+2*a(2))+exp(a(1)-2*a(2))+exp(-a(1))
   minimize(func)
cvx_end
    