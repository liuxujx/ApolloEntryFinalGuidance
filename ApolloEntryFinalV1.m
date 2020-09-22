function ApolloEntryFinalV1()
%APOLLOENTRYFINAL Apollo��������ĩ���Ƶ�
% �ο�Apollo-derived Mars precision lander guidance��֣��ԣ��ʿ���Ľ��б�д
% �����˺����Ƶ������ƺ������Ƶ���Ч����̫��
close all;
R0=6378136.3; %m
g0=9.80665; %m/s^2
GM=g0*R0^2;
% Vscale=sqrt(R0*g0); %7908.4m/s
% tscale=sqrt(R0/g0); %806.329s
W0=7.2921151e-5; % 7.2921151e-5;
J2=1.08262668e-3; % 1.08262668e-3;
rho0=1.225;
hs=7200;
d2r=pi/180;
r2d=180/pi;
mass=6800;
CL  = 0.3892;%��������ϵ��
CD  = 1.3479;%��������ϵ��
Sr   = 12.88 ; %�������ο�������m2   12.88
% ��ʼ״̬
r0 = R0+80e3; %��ʼ������ľ�
theta0 = 0*d2r; %��ʼ���ȣ����ȣ�%����ʱ����һ��ҲҪ����
phi0 = 0*d2r; %��ʼγ��
V0   = 7800; %10660; %m/s��ʼ�����ٶ� 7450
gamma0 = -1.8*d2r; %-5.8*d2r��ʼ�����
psi0 = 90*d2r; %��ʼ����ǣ�����˳ʱ��Ϊ�� 60
s0 = 0; %��ʼ����

%��ȡ�ο��켣
RefTraj=load('RefTrajV1.txt', '-ascii');%�ɺ��� RefTraj.m ���ɣ�Ϊ������
time_ref=RefTraj(:,1);
state_ref=RefTraj(:,2:8);
ctrl_ref(:,1)=RefTraj(:,9); % bank rad
ctrl_ref(:,2)=RefTraj(:,10); % �������������ٶ�L/m
ctrl_ref(:,3)=RefTraj(:,11); % �������������ٶ�D/m
ctrl_ref(:,4)=RefTraj(:,12); % ������ L/D*cos(sigma)
gammaF_ref=state_ref(end,5); % �ο��켣�����һ��gamma
pp=length(time_ref);
dt_ref=-1; %�Ͳο��켣�еĻ��ֲ�������һ��
lambda(pp,:)=[1, 0, 0, -cot(gammaF_ref), 0];
while pp>1
    %������������ٵİ��淽��
    lambda(pp-1,:)=AdjointRK4(@AdjointDyna,time_ref(pp,1),lambda(pp,:),state_ref(pp,:),ctrl_ref(pp,:),dt_ref);
    pp=pp-1;
end
%�������
ni=1; nf=1; % ���ؿ������
for n=ni:1:nf
p=1;
state(p,:)=[r0+unifrnd(-150,150), theta0+normrnd(0,0.1/3)*d2r, phi0+normrnd(0,0.1/3)*d2r,...
    V0+normrnd(0,20/3), gamma0+normrnd(0,0.1/3)*d2r, psi0+unifrnd(-0.05,0.05)*d2r, s0];
dt=1;
time(p,1)=0;
h_thres=10e3;
while 1 % 
%     r=state(p,1); V=state(p,4);%������
%     L=rho0*exp(-(r-1)*R0/hs)*(Vscale*V)^2*Sr*CL/(2*mass*g0);%�������������ٶ�
%     D=rho0*exp(-(r-1)*R0/hs)*(Vscale*V)^2*Sr*CD/(2*mass*g0);%�������������ٶ�
%     ctrl(p,1)=0*d2r;
    [ctrl(p,1),~,~,~,~,Zs(p,1),Y(p,1)]=ApolloGuidanceFinal(p,state(p,:),state_ref,ctrl_ref,lambda,theta0,phi0);
%     ctrl(p,2)=L;
%     ctrl(p,3)=D;
    state_temp=RK4(@EntryDyna3,time(p,1),state(p,:),ctrl(p,1),dt);
    h=state_temp(1)-R0; %�ȳ��Ի���һ�����ڣ��ﵽ��ɡ�߶Ⱥ�ֹͣ
    if h<=h_thres
%     V=state_temp(4);%�ȳ��Ի���һ�����ڣ��ﵽ��ɡ�ٶȺ�ֹͣ
%     if V<=state_ref(end,4)
        break;
        %break����ֱ�������ò�ѭ�� continue����ֱ�ӽ���ò�ѭ������һ�ε���
        %return����ֱ���˳������������ ��ŵĹ�ϵΪreturn>break>continue
    end
    state(p+1,:)=state_temp;
    time(p+1,1)=time(p,1)+dt;
    p=p+1;
%     pause(0.05)
end
longitude_f(n)=state(end,2)*r2d;
latitude_f(n)=state(end,3)*r2d;
end
if nf==1
    %���ӻ�
    rad=state(:,1); rad_ref=state_ref(:,1);
    speed=state(:,4);
    alt=(rad-R0)/1000; alt_ref=(rad_ref-R0)/1000;
    lon=state(:,2); lon_ref=state_ref(:,2); 
    lat=state(:,3); lat_ref=state_ref(:,3);
    figure; hold on;
    plot(time,alt,'k','linewidth',1.5)
    plot(time_ref,alt_ref,'b--','linewidth',1.5)
    xlabel('ʱ�� (s)'); ylabel('�߶� (km)')
    legend('�Ƶ�','�ο�');

    figure; hold on;
    plot3(lon*r2d,lat*r2d,alt,'k','linewidth',1.5);
    plot3(lon_ref*r2d,lat_ref*r2d,alt_ref,'b--','linewidth',1.5);
    xlabel('���� (deg)'); ylabel('γ�� (deg)'); zlabel('�߶� (km)');
    title('Three-dimensional flight path');
    legend('�Ƶ�','�ο�');

    figure; hold on;
    plot(time,ctrl*r2d,'k','linewidth',1.5);
    xlabel('ʱ�� (s)'); ylabel('���� (deg)')

    figure; hold on;
    plot(speed,Zs,'k','linewidth',1.5);
    plot(speed,[Y,-Y],'r','linewidth',1.5);
    xlabel('�ٶ� (m/s)'); ylabel('��� (m)')

    figure; hold on;
    % earthradius = almanac('earth','radius','km');
    [lat5,lon5] = scircle1(lat_ref(end),lon_ref(end),5e3,[],R0,'radians'); %Ĭ�ϽǶȣ��ĳɻ���
    [lat10,lon10] = scircle1(lat_ref(end),lon_ref(end),1e4,[],R0,'radians'); %Ĭ�ϽǶȣ��ĳɻ���
    plot(lon5*r2d,lat5*r2d,'r','linewidth',1.5);
    plot(lon10*r2d,lat10*r2d,'b','linewidth',1.5);
    plot(lon_ref(end)*r2d,lat_ref(end)*r2d,'ro','linewidth',1.5);
    plot(lon(end)*r2d,lat(end)*r2d,'kx','linewidth',1.5);
    xlabel('���� (deg)'); ylabel('γ�� (deg)');
    legend('5km','10km','Ŀ�����','ʵ�����');

    % outp=[time state ctrl];
    % save 'ActTraj.txt' outp -ascii;
else
    [lat5,lon5] = scircle1(state_ref(end,3),state_ref(end,2),5e3,[],R0,'radians'); %Ĭ�ϽǶȣ��ĳɻ���
    [lat10,lon10] = scircle1(state_ref(end,3),state_ref(end,2),1e4,[],R0,'radians'); %Ĭ�ϽǶȣ��ĳɻ���
    figure; hold on;
    plot(lon5*r2d,lat5*r2d,'r','linewidth',1.5);
    plot(lon10*r2d,lat10*r2d,'b','linewidth',1.5);
    plot(state_ref(end,2)*r2d,state_ref(end,3)*r2d,'ro','linewidth',1.5);
    plot(longitude_f,latitude_f,'kx','linewidth',1.5)
end

end

function [ctrl,d_u,F1,F2,F3,Zs,Y]=ApolloGuidanceFinal(p,state,state_ref,ctrl_ref,lambda,theta0,phi0)
persistent sign_ctrl

R0=6378136.3; %m
g0=9.80665; %m/s^2
GM=g0*R0^2;
Vscale=sqrt(R0*g0); %7908.4m/s
tscale=sqrt(R0/g0); %806.329s
W0=7.2921151e-5; % 7.2921151e-5;
J2=1.08262668e-3; % 1.08262668e-3;
rho0=1.225;
hs=7200;
d2r=pi/180;
r2d=180/pi;
mass=6800;
CL  = 0.3892;%��������ϵ��
CD  = 1.3479;%��������ϵ��
Sr   = 12.88 ; %�������ο�������m2   12.88
%Apolloĩ�Ƶ�
%������״̬ state=[r, theta, phi, V, gamma, psi, s]
r=state(1); theta=state(2); phi=state(3);
V=state(4); gamma=state(5); psi=state(6); s=state(7);
L=rho0*exp(-(r-R0)/hs)*V^2*Sr*CL/2;
D=rho0*exp(-(r-R0)/hs)*V^2*Sr*CD/2;
r_interp=interp1(state_ref(:,4),state_ref(:,1),V,'linear','extrap');
gamma_interp=interp1(state_ref(:,4),state_ref(:,5),V,'linear','extrap');
s_interp=interp1(state_ref(:,4),state_ref(:,7),V,'linear','extrap');
ctrl_interp=interp1(state_ref(:,4),ctrl_ref(:,1),V,'linear','extrap');
L_interp=interp1(state_ref(:,4),ctrl_ref(:,2),V,'linear','extrap');
D_interp=interp1(state_ref(:,4),ctrl_ref(:,3),V,'linear','extrap');
u_interp=interp1(state_ref(:,4),ctrl_ref(:,4),V,'linear','extrap');
%lambda=[lambda_s, lambda_V, lambda_gamma, lambda_r, lambda_u]
lambda_r_interp=interp1(state_ref(:,4),lambda(:,4),V,'linear','extrap');
lambda_gamma_interp=interp1(state_ref(:,4),lambda(:,3),V,'linear','extrap');
lambda_u_interp=interp1(state_ref(:,4),lambda(:,5),V,'linear','extrap');
F1=lambda_r_interp/(-rho0*exp(-(r_interp-R0)/hs)*V^2*Sr*CD/(2*hs)); % *mass
F2=lambda_gamma_interp/(V*cos(gamma_interp));
F3=lambda_u_interp;
[F1 F2 F3];
K=1;%����ϵ�� over control,��������K=5 ���ǽ���K=2
% ״̬ƫ��
d_s=s-s_interp;
d_rdot=V*sin(gamma)-V*sin(gamma_interp);
d_D=D-D_interp;
s_bias=0e3; % s_bias ������̬����ϵͳ�����Ƿ�ת���µ��ݳ����
d_u=K/F3*(d_s+F2*d_rdot+F1*d_D+s_bias);
% u=L_interp/D_interp*cos(ctrl_interp)+d_u; %u=L/D*cos(sigma)
u=u_interp+d_u;
if u*D/L>cosd(10)
    ctrl=10*d2r;
    'Waring: abs(sigma)<10 '
elseif u*D/L<cosd(170)
    ctrl=170*d2r;
    'Waring: abs(sigma)>170 '
else
    ctrl=acos(u*D/L); %[0 pi]
end

% �����Ƶ�
% ��̺��ݳ̼��� Advances in spacecraft atmospheric entry guidance (Page43)
thetaF=state_ref(end,2);
phiF=state_ref(end,3);
deltae0=acos(sin(phi0)*sin(phiF)+cos(phi0)*cos(phiF)*cos(theta0-thetaF)); % (2.40) d12
deltae=acos(sin(phi0)*sin(phi)+cos(phi0)*cos(phi)*cos(theta0-theta)); % (2.40) d13
c0=sign(thetaF-theta0)*acos((sin(phiF)-sin(phi0)*cos(deltae0))/cos(phi0)/sin(deltae0)); % (2.43)
sig0=pi/2-c0; %��ʼ���䷽λ�ǣ��Ӷ�����ʱ�� Psi12
cs=sign(theta-theta0)*acos((sin(phi)-sin(phi0)*cos(deltae))/cos(phi0)/sin(deltae)); % (2.43)
sigs=pi/2-cs; %˲ʱ���䷽λ�� Psi13
deltaZs=asin(sin(deltae)*sin(sig0-sigs)); % (2.45)
deltaRs=acos(cos(deltae)/cos(deltaZs)); % (2.46)

Zs=R0*deltaZs;  %�����̴�С
Rs=R0*deltaRs; %�����ݳ̴�С
% ©���߽�
V1=150; K1=600; Y1=1e3;
if V>3600
    Y=15e3;
else
    Y=(V-V1)^2/K1+Y1;
end
% ���ŷ�ת
if p==1
    sign_ctrl=1;
end  
if abs(Zs)>=Y
%     ctrl=-abs(ctrl)*sign_ctrl;
    ctrl=-ctrl*sign_ctrl;
else
    ctrl=ctrl*sign_ctrl;
end

end

function outp=AdjointDyna(t,inp1,inp2,ctrl)
R0=6378136.3; %m
g0=9.80665; %m/s^2
GM=g0*R0^2;
% Vscale=sqrt(R0*g0); %7908.4m/s
% tscale=sqrt(R0/g0); %806.329s
W0=7.2921151e-5; % 7.2921151e-5
J2=1.08262668e-3; % 1.08262668e-3;
rho0=1.225;
hs=7200;
d2r=pi/180;
r2d=180/pi;
mass=6800;
CL  = 0.3892;%��������ϵ��
CD  = 1.3479;%��������ϵ��
Sr   = 12.88 ; %�������ο�������m2   12.88

%AdjointDyna ���淽��
% inp1=[lambda_s, lambda_V, lambda_gamma, lambda_r, lambda_u]
% inp2=[r_ref, theta_ref, phi_ref, V_ref, gamma_ref, psi_ref, s_ref]
% ctrl=sigma_ref; �ο��켣����ʱ���ֵ�õ���������״̬������
s=inp2(7); V=inp2(4); gamma=inp2(5); r=inp2(1);
sigma=ctrl(1); u=ctrl(4); % L/D*cos(sigma)
A=[0, cos(gamma), -V*sin(gamma), 0; %g=GM/r^2������g��r�ĵ���
    0, -rho0*exp(-(r-R0)/hs)*V*Sr*CD/mass, -GM/r^2*cos(gamma),...
    rho0*exp(-(r-R0)/hs)*V^2*Sr*CD/(2*mass*hs)+2*GM/r^3*sin(gamma);
    0, rho0*exp(-(r-R0)/hs)*Sr*CD*u/(2*mass)+(1/r+GM/(r^2*V^2))*cos(gamma),...
    -(V/r-GM/(r^2*V))*sin(gamma),...
    -rho0*exp(-(r-R0)/hs)*V*Sr*CD*u/(2*mass*hs)+(-V/r^2+2*GM/(r^3*V))*cos(gamma);
    0, sin(gamma), V*cos(gamma), 0];
% A=[0, cos(gamma), -V*sin(gamma), 0; %g=GM/r^2��������g��r�ĵ���
%     0, -rho0*exp(-(r-R0)/hs)*V*Sr*CD/mass, -g0*cos(gamma),...
%     rho0*exp(-(r-R0)/hs)*V^2*Sr*CD/(2*mass*hs);
%     0, rho0*exp(-(r-R0)/hs)*Sr*CL*cos(sigma)/(2*mass)+(1/r+g0/V^2)*cos(gamma),...
%     -(V/r-g0/V)*sin(gamma),...
%     -rho0*exp(-(r-R0)/hs)*V*Sr*CL*cos(sigma)/(2*mass*hs)+(-V/r^2)*cos(gamma);
%     0, sin(gamma), V*cos(gamma), 0];
AA=[[-A.' zeros(4,1)]; [0, 0, rho0*exp(-(r-R0)/hs)*V*Sr*CD/(2*mass), 0, 0]];%��lambda�������� u=L/D*cos(sigma)
% AA=[[-A.' zeros(4,1)]; [0, 0, rho0*exp(-(r-R0)/hs)*V*Sr*CL/(2*mass), 0, 0]];%��lambda�������� u=cos(sigma)
outp=(AA*inp1.').';
end

function outp=AdjointRK4(func,t,inp1,inp2,ctrl,dt)
k1=func(t,inp1,inp2,ctrl);
k2=func(t+dt/2,inp1+dt/2*k1,inp2,ctrl);
k3=func(t+dt/2,inp1+dt/2*k2,inp2,ctrl);
k4=func(t+dt,inp1+dt*k3,inp2,ctrl);
outp=inp1+ dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function outp=EntryDyna3(t,inp,ctrl)
R0=6378136.3; %m
g0=9.80665; %m/s^2
GM=g0*R0^2;
Vscale=sqrt(R0*g0); %7908.4m/s
tscale=sqrt(R0/g0); %806.329s
W0=7.2921151e-5; % 7.2921151e-5;
J2=1.08262668e-3; % 1.08262668e-3;
rho0=1.225;
hs=7200;
d2r=pi/180;
r2d=180/pi;
mass=6800;
CL  = 0.3892;%��������ϵ��
CD  = 1.3479;%��������ϵ��
Sr   = 12.88 ; %�������ο�������m2   12.88

%�����ٲ���
r     = inp(1); %���ľ�
theta = inp(2); %����
phi   = inp(3); %γ��
V     = inp(4); %�ٶ�
gamma = inp(5); %������� ����·����
psi   = inp(6); %����ƫ�� ������λ��
s     = inp(7); %����
sigma = ctrl; %����
%��ر���
h=r-R0;
L=rho0*exp(-h/hs)*V^2*Sr*CL/2;
D=rho0*exp(-h/hs)*V^2*Sr*CD/2;
%�����ٶ���ѧģ��
rdot = V*sin(gamma);
thetadot = V*cos(gamma)*sin(psi)/(r*cos(phi));
phidot = V*cos(gamma)*cos(psi)/r;
Vdot = -D/mass-GM/r^2*sin(gamma);
gammadot = L*cos(sigma)/(mass*V)-(GM/(r^2*V)-V/r)*cos(gamma);
psidot = L*sin(sigma)/(mass*V*cos(gamma))+V*cos(gamma)*sin(psi)*tan(phi)/r;
sdot = V*cos(gamma);
outp=[rdot thetadot phidot Vdot gammadot psidot sdot];
end

function outp=RK4(func,t,inp,ctrl,dt)
k1=func(t,inp,ctrl);
k2=func(t+dt/2,inp+dt/2*k1,ctrl);
k3=func(t+dt/2,inp+dt/2*k2,ctrl);
k4=func(t+dt,inp+dt*k3,ctrl);
outp=inp+ dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

