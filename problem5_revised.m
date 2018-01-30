%linear algebra term problem 5
%A=[0.2 0.5;0.3 0.7;1 2;1 3;0.9 1.4;0.4 0.3;1.1 1.2;0.8 2.2];
%B=[-1.1 2;-0.3 0.2;-0.7 0.9;-2 0.4;-1.3 1.4;-2.2 2.4;-0.7 1.7;-1.4 .6];
A=randn(8,2);
B=randn(8,2);
a=[0.5 0.5 0.5];
A=[A -ones(8,1)];
B=[B -ones(8,1)];
u=1.5;
tbarrier=1.2;
maxiter=20;
stop=0.01;
alpha=0.1;
beta=0.5;
value=[];
for j=1:8
    value=[value a*A(j,:)' a*B(j,:)'];
end
value=[value(1:2:15);value(2:2:16)];
if -min(value(1,:))>max(value(2,:))
    S=-min(value(1,:))+0.1;
else
    S=max(value(2,:))+0.1;
end
for i=1:maxiter
    %gradient method
    sum1=0;
    sum2=0;
    sum3=0;
    for j=1:8
        sum1=sum1-A(j,1)/(a*A(j,:)'+S)+B(j,1)/(-a*B(j,:)'+S);
        sum2=sum2-A(j,2)/(a*A(j,:)'+S)+B(j,2)/(-a*B(j,:)'+S);
        sum3=sum3+1/(a*A(j,:)'+S)-1/(-a*B(j,:)'+S);
    end
    gradient=-[sum1;sum2;sum3];    %gradient step
    while 1      %stop criterion of gradient descent method
        %backtracking line search
        func=tbarrier*S;
        for j=1:8
            func=func-log(a*A(j,:)'+S)-log(-a*B(j,:)'+S);
        end
        tline=1;
        atemp=a+tline*gradient';               %A+t*deltaA
        value=[];
        functemp=tbarrier*S;
        for j=1:8
            functemp=functemp-log(atemp*A(j,:)'+S)-log(-atemp*B(j,:)'+S);
        end
        while functemp>func+alpha*tline*-gradient'*gradient
            tline=beta*tline;
            atemp=a+tline*gradient';
            functemp=tbarrier*S;
            for j=1:8
                functemp=functemp-log(atemp*A(j,:)'+S)-log(-atemp*B(j,:)'+S);
            end
        end
        a=a+tline*gradient';
        sum1=0;
        sum2=0;
        sum3=0;
        for j=1:8
            sum1=sum1-A(j,1)/(a*A(j,:)'+S)+B(j,1)/(-a*B(j,:)'+S);
            sum2=sum2-A(j,2)/(a*A(j,:)'+S)+B(j,2)/(-a*B(j,:)'+S);
            sum3=sum3+1/(a*A(j,:)'+S)-1/(-a*B(j,:)'+S);
        end
        gradient=-[sum1;sum2;sum3];
        if norm(gradient)<0.01
            break;
        end
        value=[];
        for j=1:8
            value=[value a*A(j,:)' a*B(j,:)'];
        end
        value=[value(1:2:15);value(2:2:16)];
        if -min(value(1,:))>max(value(2,:))
            S=-min(value(1,:))+0.1;
        else
            S=max(value(2,:))+0.1;
        end
        vals=tbarrier*S;
        for j=1:8
            vals=vals-log(a*A(j,:)'+S)-log(-a*B(j,:)'+S);
        end
        if vals>func
            break;
        end;
        if S<0
            break;
        end
    end
        if (17/tbarrier)<stop
            break;
        end
    tbarrier=tbarrier*u;
    value=[];
    for j=1:8
        value=[value a*A(j,:)' a*B(j,:)'];
    end
    value=[value(1:2:15);value(2:2:16)];
    if -min(value(1,:))>max(value(2,:))
        S=-min(value(1,:))+0.1;
    else
        S=max(value(2,:))+0.1;
    end
    if S<0
        break;
    end
end
 plot(B(:,1),B(:,2),'o');
 hold;
 plot(A(:,1),A(:,2),'+');
 x=-3:3;
 plot(x,a(1)*x-a(3));
 hold;