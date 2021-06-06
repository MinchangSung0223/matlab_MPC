close all;
g  = 9.8; % gravity
zc = 0.8; % Center of Mass Height
Ts = 5e-3; % Sampling Time
EndTime = 10;
ts = linspace(0,EndTime,EndTime/Ts)



A = [0 1 0;...
     0 0 1;...
     0 0 0]
B = [0;0;1];
C= [1 0 -zc/g; 1 0 0; 0 1 0];
D = [0;0;0];
lip_x = ss(A,B,C,D);
lip_y = ss(A,B,C,D);
lip_x_d = c2d(lip_x,Ts);
lip_y_d = c2d(lip_y,Ts);
lip_x_d=setmpcsignals(lip_x_d,"MeasuredOutputs",1,"UnmeasuredOutputs",[2 3]);
lip_y_d=setmpcsignals(lip_y_d,"MeasuredOutputs",1,"UnmeasuredOutputs",[2 3]);
xm = [0;0;0];


n1 = size(xm,1);
q = 1;
Np = 300;
Nc = 100;


Am = lip_x_d.A;
Bm = lip_x_d.B;
Cm = lip_x_d.C;
Cm = Cm(1,:);
Dm = lip_x_d.D;
y = Cm*xm;

xm_old = xm;
A = [Am zeros(n1,q);Cm*Am eye(q)];
B = [Bm;Cm*Bm];
C = [zeros(q,n1) eye(q,q)];
x = [xm-xm_old;y]

u=zeros(size(ts))
u(200:end)=1;
k = 1;
data=[]


 F = []
 Phi= []
 for i =1:1:Np
     F =[F;C*A^i];
     temp =[]
     for j = 1:1:Nc
        if i-j<0
             temp = [temp,zeros(size(C*B))];
        else
             temp = [temp,C*A^(i-j)*B];
        end
     end
     Phi =[Phi;temp];
 end
rw = 0.01;
R_bar = rw*eye(Nc);
r = u;
R_bar_s = ones(size(r,1),Np*q);
Phi_Phi = Phi'*Phi;
Phi_R = Phi'*R_bar_s';
Phi_F = Phi'*F;
u=0;
Kob =1;
xhat = [0;0;0];
for t = ts
    DeltaU = inv(Phi_Phi+rw*eye(Nc,Nc))*(Phi_R*r(k) -Phi_F*x);
    deltau = DeltaU(1,1);
    u = u+deltau;
    xm_old = xm;
    xm = Am*xm+Bm*u;
    y = Cm*xm;
    xhat = Am*xhat+Bm*u+Kob*(y(1)-xhat(1))
    x = [xm-xm_old;y(1)];

    
     data(1,k) = t;
     data(2,k) = xm(1);
     data(3,k) = xm(2);
     data(4,k) = xm(3);
     data(5,k) = xhat(1);
     data(6,k) = xhat(2);
     data(7,k) = xhat(3);
     
     
     data(8,k) = u;
     k = k+1;
end


figure;
subplot(4,1,1)
plot(data(1,:),data(2,:))
hold on;
plot(data(1,:),data(5,:),'r:')

subplot(4,1,2)
plot(data(1,:),data(3,:))
hold on;
plot(data(1,:),data(6,:),'r:')
subplot(4,1,3)
plot(data(1,:),data(4,:))
hold on;
plot(data(1,:),data(7,:),'r:')
subplot(4,1,4)
plot(data(1,:),data(8,:))
