close all;
g  = 9.8; % gravity
zc = 0.8; % Center of Mass Height
Ts = 5e-3; % Sampling Time
EndTime = 10;
sim_time = linspace(0,EndTime,EndTime/Ts)



Ac = [0 1 0;...
     0 0 1;...
     0 0 0]
Bc = [0;0;1];
Cc= [1 0 -zc/g];
Dc = [0];

lip_x = ss(Ac,Bc,Cc,Dc);
lip_y = ss(Ac,Bc,Cc,Dc);
lip_x_d = c2d(lip_x,Ts);
lip_y_d = c2d(lip_y,Ts);
lip_x_d=setmpcsignals(lip_x_d,"MeasuredOutputs",1);
lip_y_d=setmpcsignals(lip_y_d,"MeasuredOutputs",1);

Ap = lip_x_d.A;
Bp = lip_x_d.B;
Cp = lip_x_d.C;
Dp = lip_x_d.D;
Nc = 50;
Np = 300;
Nl = 200;
rw =0.00005;
BarR = rw*eye(Nc);
[Phi_Phi, Phi_F, Phi_R, F, BarRs, Phi,Psi, A_e, B_e,C_e]= mpcgain(Ap, Bp,zeros(3,1), Cp, Nc, Np,Nl);

k= 1;
xm = [0;0;0];
old_xm = [0;0;0];
data = []
y = Cp*xm;
x = [xm;y]
r = zeros(size(sim_time,2)+Nl,1)
r(size(sim_time,2)/7:size(sim_time,2)*2/7) = 1;
r(size(sim_time,2)*2/7:size(sim_time,2)*3/7) = -1;
r(size(sim_time,2)*3/7:size(sim_time,2)*4/7) = 1;
r(size(sim_time,2)*4/7:size(sim_time,2)*5/7) = -1;


r_diff = []
for i = 2:1:length(r)
    r_diff = [r_diff,r(i)-r(i-1)];
end
r_diff(end+1) = r_diff(end);

u =0;
Q = eye(Np,Np);
for i = 1:1:Np
    Q(i,i) = 1/log10(1+i);
end
for t = sim_time
    DeltaU  = -inv(Phi'*Q*Phi+BarR)*(Phi'*Q*Psi*r_diff(k+1:k+1+Nl-1)' +Phi'*Q*F*x);
    deltau = DeltaU(1,1);
    u = u+deltau;
    old_xm = xm;
    xm = Ap*xm+Bp*u;
    y = Cp*xm;
    x = [xm-old_xm;y];
    
    data(k,1) = t;
    data(k,2) = y(1);
    data(k,3) = xm(2);
    data(k,4) = xm(3);
    data(k,5) = u;
    data(k,6) = r(k);    
    k = k+1;
end

plot_data(data);
