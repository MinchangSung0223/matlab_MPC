function [Phi_Phi, Phi_F, Phi_R, F, BarRs, Phi,Psi, A_e, B_e,C_e,B_p]= mpcgain(Ap, Bp,Bp2, Cp, Nc, Np,Nl)
[m1,n1] = size(Cp);
[n1,n_in] = size(Bp);
[n12,n_in2] = size(Bp2);

A_e = eye(n1+m1,n1+m1);
A_e(1:n1,1:n1) = Ap;
A_e(n1+1:n1+m1,1:n1) = Cp * Ap;
B_e = zeros(n1+m1,n_in);
B_p = zeros(n12+m1,n_in2);
B_p(end,:) = -1;
B_e(1:n1,:) = Bp;

B_e(n1+1:n1+m1,:) = Cp * Bp;

C_e = zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1) = eye(m1,m1);

n = n1 + m1;
h(1,:) = C_e;
F(1,:) = C_e * A_e;
for kk=2:Np
h(kk,:) = h(kk-1,:) * A_e;
F(kk,:) = F(kk-1,:) * A_e;
end
v = h * B_e;
v2 = h* B_p;
Phi = zeros(Np,Nc); %declare the dimension of Phi
Psi =  zeros(Np,Nl); %declare the dimension of Psi
Phi(:,1) = v; % first column of Phi
Psi(:,1) = v2; % first column of Psi

for i=2:Nc
Phi(:,i) = [zeros(i-1,1); v(1:Np-i+1,1)]; %Toeplitz matrix
end
for i=2:Nl
Psi(:,i) = [zeros(i-1,1); v2(1:Np-i+1,1)]; %Toeplitz matrix
end

BarRs = ones(Np,1);
Phi_Phi = Phi' * Phi;
Phi_F = Phi' * F;
Phi_R = Phi' * BarRs;
end
