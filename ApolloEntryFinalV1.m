function ApolloEntryFinalV1()
%APOLLOENTRYFINAL Apollo大气再入末段制导
% 参考Apollo-derived Mars precision lander guidance和郑艺裕博士论文进行编写
% 考虑了横向制导，但似乎横向制导的效果不太好
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
CL  = 0.3892;%气动升力系数
CD  = 1.3479;%气动阻力系数
Sr   = 12.88 ; %航天器参考横截面积m2   12.88
% 初始状态
r0 = R0+80e3; %初始再入地心距
theta0 = 0*d2r; %初始经度（弧度）%更改时下面一处也要更改
phi0 = 0*d2r; %初始纬度
V0   = 7800; %10660; %m/s初始再入速度 7450
gamma0 = -1.8*d2r; %-5.8*d2r初始再入角
psi0 = 90*d2r; %初始航向角，正北顺时针为正 60
s0 = 0; %初始航程

%读取参考轨迹
RefTraj=load('RefTrajV1.txt', '-ascii');%由函数 RefTraj.m 生成，为有量纲
time_ref=RefTraj(:,1);
state_ref=RefTraj(:,2:8);
ctrl_ref(:,1)=RefTraj(:,9); % bank rad
ctrl_ref(:,2)=RefTraj(:,10); % 有量纲升力加速度L/m
ctrl_ref(:,3)=RefTraj(:,11); % 有量纲阻力加速度D/m
ctrl_ref(:,4)=RefTraj(:,12); % 控制量 L/D*cos(sigma)
gammaF_ref=state_ref(end,5); % 参考轨迹的最后一个gamma
pp=length(time_ref);
dt_ref=-1; %和参考轨迹中的积分步长保持一致
lambda(pp,:)=[1, 0, 0, -cot(gammaF_ref), 0];
while pp>1
    %反向积分有量纲的伴随方程
    lambda(pp-1,:)=AdjointRK4(@AdjointDyna,time_ref(pp,1),lambda(pp,:),state_ref(pp,:),ctrl_ref(pp,:),dt_ref);
    pp=pp-1;
end
%仿真参数
ni=1; nf=1; % 蒙特卡洛仿真
for n=ni:1:nf
p=1;
state(p,:)=[r0+unifrnd(-150,150), theta0+normrnd(0,0.1/3)*d2r, phi0+normrnd(0,0.1/3)*d2r,...
    V0+normrnd(0,20/3), gamma0+normrnd(0,0.1/3)*d2r, psi0+unifrnd(-0.05,0.05)*d2r, s0];
dt=1;
time(p,1)=0;
h_thres=10e3;
while 1 % 
%     r=state(p,1); V=state(p,4);%有量纲
%     L=rho0*exp(-(r-1)*R0/hs)*(Vscale*V)^2*Sr*CL/(2*mass*g0);%无量纲升力加速度
%     D=rho0*exp(-(r-1)*R0/hs)*(Vscale*V)^2*Sr*CD/(2*mass*g0);%无量纲阻力加速度
%     ctrl(p,1)=0*d2r;
    [ctrl(p,1),~,~,~,~,Zs(p,1),Y(p,1)]=ApolloGuidanceFinal(p,state(p,:),state_ref,ctrl_ref,lambda,theta0,phi0);
%     ctrl(p,2)=L;
%     ctrl(p,3)=D;
    state_temp=RK4(@EntryDyna3,time(p,1),state(p,:),ctrl(p,1),dt);
    h=state_temp(1)-R0; %先尝试积分一个周期，达到开伞高度后停止
    if h<=h_thres
%     V=state_temp(4);%先尝试积分一个周期，达到开伞速度后停止
%     if V<=state_ref(end,4)
        break;
        %break就是直接跳出该层循环 continue就是直接进入该层循环的下一次迭代
        %return就是直接退出程序或函数返回 大概的关系为return>break>continue
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
    %可视化
    rad=state(:,1); rad_ref=state_ref(:,1);
    speed=state(:,4);
    alt=(rad-R0)/1000; alt_ref=(rad_ref-R0)/1000;
    lon=state(:,2); lon_ref=state_ref(:,2); 
    lat=state(:,3); lat_ref=state_ref(:,3);
    figure; hold on;
    plot(time,alt,'k','linewidth',1.5)
    plot(time_ref,alt_ref,'b--','linewidth',1.5)
    xlabel('时间 (s)'); ylabel('高度 (km)')
    legend('制导','参考');

    figure; hold on;
    plot3(lon*r2d,lat*r2d,alt,'k','linewidth',1.5);
    plot3(lon_ref*r2d,lat_ref*r2d,alt_ref,'b--','linewidth',1.5);
    xlabel('经度 (deg)'); ylabel('纬度 (deg)'); zlabel('高度 (km)');
    title('Three-dimensional flight path');
    legend('制导','参考');

    figure; hold on;
    plot(time,ctrl*r2d,'k','linewidth',1.5);
    xlabel('时间 (s)'); ylabel('倾侧角 (deg)')

    figure; hold on;
    plot(speed,Zs,'k','linewidth',1.5);
    plot(speed,[Y,-Y],'r','linewidth',1.5);
    xlabel('速度 (m/s)'); ylabel('横程 (m)')

    figure; hold on;
    % earthradius = almanac('earth','radius','km');
    [lat5,lon5] = scircle1(lat_ref(end),lon_ref(end),5e3,[],R0,'radians'); %默认角度，改成弧度
    [lat10,lon10] = scircle1(lat_ref(end),lon_ref(end),1e4,[],R0,'radians'); %默认角度，改成弧度
    plot(lon5*r2d,lat5*r2d,'r','linewidth',1.5);
    plot(lon10*r2d,lat10*r2d,'b','linewidth',1.5);
    plot(lon_ref(end)*r2d,lat_ref(end)*r2d,'ro','linewidth',1.5);
    plot(lon(end)*r2d,lat(end)*r2d,'kx','linewidth',1.5);
    xlabel('经度 (deg)'); ylabel('纬度 (deg)');
    legend('5km','10km','目标落点','实际落点');

    % outp=[time state ctrl];
    % save 'ActTraj.txt' outp -ascii;
else
    [lat5,lon5] = scircle1(state_ref(end,3),state_ref(end,2),5e3,[],R0,'radians'); %默认角度，改成弧度
    [lat10,lon10] = scircle1(state_ref(end,3),state_ref(end,2),1e4,[],R0,'radians'); %默认角度，改成弧度
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
CL  = 0.3892;%气动升力系数
CD  = 1.3479;%气动阻力系数
Sr   = 12.88 ; %航天器参考横截面积m2   12.88
%Apollo末制导
%有量纲状态 state=[r, theta, phi, V, gamma, psi, s]
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
K=1;%过控系数 over control,地球再入K=5 火星进入K=2
% 状态偏差
d_s=s-s_interp;
d_rdot=V*sin(gamma)-V*sin(gamma_interp);
d_D=D-D_interp;
s_bias=0e3; % s_bias 补偿姿态控制系统因倾侧角翻转导致的纵程误差
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

% 横向制导
% 横程和纵程计算 Advances in spacecraft atmospheric entry guidance (Page43)
thetaF=state_ref(end,2);
phiF=state_ref(end,3);
deltae0=acos(sin(phi0)*sin(phiF)+cos(phi0)*cos(phiF)*cos(theta0-thetaF)); % (2.40) d12
deltae=acos(sin(phi0)*sin(phi)+cos(phi0)*cos(phi)*cos(theta0-theta)); % (2.40) d13
c0=sign(thetaF-theta0)*acos((sin(phiF)-sin(phi0)*cos(deltae0))/cos(phi0)/sin(deltae0)); % (2.43)
sig0=pi/2-c0; %初始发射方位角，从东向逆时针 Psi12
cs=sign(theta-theta0)*acos((sin(phi)-sin(phi0)*cos(deltae))/cos(phi0)/sin(deltae)); % (2.43)
sigs=pi/2-cs; %瞬时发射方位角 Psi13
deltaZs=asin(sin(deltae)*sin(sig0-sigs)); % (2.45)
deltaRs=acos(cos(deltae)/cos(deltaZs)); % (2.46)

Zs=R0*deltaZs;  %计算横程大小
Rs=R0*deltaRs; %计算纵程大小
% 漏斗边界
V1=150; K1=600; Y1=1e3;
if V>3600
    Y=15e3;
else
    Y=(V-V1)^2/K1+Y1;
end
% 符号翻转
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
CL  = 0.3892;%气动升力系数
CD  = 1.3479;%气动阻力系数
Sr   = 12.88 ; %航天器参考横截面积m2   12.88

%AdjointDyna 伴随方程
% inp1=[lambda_s, lambda_V, lambda_gamma, lambda_r, lambda_u]
% inp2=[r_ref, theta_ref, phi_ref, V_ref, gamma_ref, psi_ref, s_ref]
% ctrl=sigma_ref; 参考轨迹按照时间插值得到的有量纲状态、控制
s=inp2(7); V=inp2(4); gamma=inp2(5); r=inp2(1);
sigma=ctrl(1); u=ctrl(4); % L/D*cos(sigma)
A=[0, cos(gamma), -V*sin(gamma), 0; %g=GM/r^2，考虑g对r的导数
    0, -rho0*exp(-(r-R0)/hs)*V*Sr*CD/mass, -GM/r^2*cos(gamma),...
    rho0*exp(-(r-R0)/hs)*V^2*Sr*CD/(2*mass*hs)+2*GM/r^3*sin(gamma);
    0, rho0*exp(-(r-R0)/hs)*Sr*CD*u/(2*mass)+(1/r+GM/(r^2*V^2))*cos(gamma),...
    -(V/r-GM/(r^2*V))*sin(gamma),...
    -rho0*exp(-(r-R0)/hs)*V*Sr*CD*u/(2*mass*hs)+(-V/r^2+2*GM/(r^3*V))*cos(gamma);
    0, sin(gamma), V*cos(gamma), 0];
% A=[0, cos(gamma), -V*sin(gamma), 0; %g=GM/r^2，不考虑g对r的导数
%     0, -rho0*exp(-(r-R0)/hs)*V*Sr*CD/mass, -g0*cos(gamma),...
%     rho0*exp(-(r-R0)/hs)*V^2*Sr*CD/(2*mass*hs);
%     0, rho0*exp(-(r-R0)/hs)*Sr*CL*cos(sigma)/(2*mass)+(1/r+g0/V^2)*cos(gamma),...
%     -(V/r-g0/V)*sin(gamma),...
%     -rho0*exp(-(r-R0)/hs)*V*Sr*CL*cos(sigma)/(2*mass*hs)+(-V/r^2)*cos(gamma);
%     0, sin(gamma), V*cos(gamma), 0];
AA=[[-A.' zeros(4,1)]; [0, 0, rho0*exp(-(r-R0)/hs)*V*Sr*CD/(2*mass), 0, 0]];%对lambda进行增广 u=L/D*cos(sigma)
% AA=[[-A.' zeros(4,1)]; [0, 0, rho0*exp(-(r-R0)/hs)*V*Sr*CL/(2*mass), 0, 0]];%对lambda进行增广 u=cos(sigma)
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
CL  = 0.3892;%气动升力系数
CD  = 1.3479;%气动阻力系数
Sr   = 12.88 ; %航天器参考横截面积m2   12.88

%无量纲参数
r     = inp(1); %地心距
theta = inp(2); %经度
phi   = inp(3); %纬度
V     = inp(4); %速度
gamma = inp(5); %航迹倾角 飞行路径角
psi   = inp(6); %航迹偏角 航迹方位角
s     = inp(7); %航程
sigma = ctrl; %倾侧角
%相关变量
h=r-R0;
L=rho0*exp(-h/hs)*V^2*Sr*CL/2;
D=rho0*exp(-h/hs)*V^2*Sr*CD/2;
%有量纲动力学模型
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

