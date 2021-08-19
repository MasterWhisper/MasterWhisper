%Lin Jian博士论文
%M*q''+Ωc*G*q+(K-Ωc^2*KΩ)*q=0   M为质量矩阵矩阵 Ωc是行星架角速度 G是陀螺矩阵 KΩ是向心矩阵 
%q是位移矩阵
clear all
clc
%% 基本参数
m=27; %模数
zs=23;
zp=17;
zr=57;
rs=27*23/1000/2;
rr=27*57/1000/2;
rp=27*17/1000/2;
rc=rs+rp;
%% 质量矩阵  18×18 M
mc=5091;
mr=900;
ms=750;
m1=525;
m2=525;
m3=525;
c=[mc,mc,5724];  %尾末最后一位是I/r^2
Mc=diag(c);
r=[mr,mr,1092];
Mr=diag(r);
s=[ms,ms,663];
Ms=diag(s);
p1=[m1,m1,420];
M1=diag(p1);
M2=M1;
M3=M2;
Mtotal=[c,r,s,p1,p1,p1];
M=diag(Mtotal);  %质量矩阵 依次为行星架 内齿圈 太阳轮 行星轮1 2 3 4 5
%% 向心矩阵  18×18 Komoge
co=[mc,mc,0];  %尾末最后一位是I/r^2
Mco=diag(co);
ro=[mr,mr,0];
Mro=diag(ro);
so=[ms,ms,0];
Mso=diag(so);
p1o=[m1,m1,0];
M1o=diag(p1o);
M2o=M1o;
M3o=M2o;
Mtotalo=[co,ro,so,p1o,p1o,p1o];
Komoge=diag(Mtotalo); %质量矩阵 依次为行星架 内齿圈 太阳轮 行星轮1 2 3 4 5
%% 支撑刚度矩阵   18×18
kcx=1e8;  %支撑刚度
kcy=1e8;
kcu=0;    %旋转刚度
kc=[kcx,kcy,kcu];
Kcb=diag(kc);
krx=1e8;  %支撑刚度
kry=1e8;
kru=1e9;    %旋转刚度
kr=[krx,kry,kru];
Krb=diag(kr);
ksx=1e8;
kxy=1e8;
ksu=0;
ks=[ksx,kxy,ksu];
Ksb=diag(ks);
kpx=1e8;
kpy=1e8;
kpu=0;
kp=[kpx,kpy,kpu];
Kp=diag(kp);
B=[kc,kr,ks,kp,kp,kp];
Kb=diag(B);
%% 陀螺矩阵 G 18×18
o=[0 0 0; 0 0 0; 0 0 0];
Gc=[0 -2*mc 0; 2*mc 0 0;0 0 0];
Gr=[0 -2*mr 0; 2*mr 0 0;0 0 0];
Gs=[0 -2*ms 0; 2*ms 0 0;0 0 0];
Gp=[0 -2*m1 0; 2*m1 0 0;0 0 0];
G=[Gc o o o o o ;
    o Gr o o o o ;
    o o Gs o o o ;
    o o o Gp o o ;
    o o o o Gp o ;
    o o o o o Gp ];
%%  内啮合时变啮合刚度
N1=17; N2=57; L=0.016; a0=20*pi/180; nv=75*pi/180;
% nv denotes angle of start point of crack
E=2.068e11; v=0.3; rintp=20; rintg=30;
%===================================
inva0=tan(a0)-a0;          a2p=pi/(2*N1)+inva0;
inva0=tan(a0)-a0;          a2g=pi/(2*N2)+inva0;
%===================================
LOA=((N2+2)^2+(N1+N2)^2-2*(N2+2)*(N1+N2)*...
    cos(acos(N2*cos(a0)/(N2+2))-a0))^0.5;
a110p=tan(acos(N1*cos(a0)/LOA))-a2p;
a110g=tan(acos(N2*cos(a0)/(N2+2)))-a2g;
%====================================
thetaD=tan(acos(N1*cos(a0)/(N1+2)))-2*pi/N1-...
    tan(acos(N1*cos(a0)/LOA));
%=======Initialization of matrix================
kb11=zeros(roundn(2*pi/N1,-3)*1000-1,1);
kb12=zeros(roundn(thetaD,-3)*1000,1);
ka11=kb11; ka21=kb11;            ka12=kb12; ka22=kb12;
ks11=kb11; ks21=kb11;            ks12=kb12; ks22=kb12;
kb21=kb11;                               kb22=kb12;
kf11=kb11; kf21=kb11;            kf12=kb12; kf22=kb12;
% First subscript 1 or 2  indecates the pinion or the gear
% Second subscript 1 indecates the first tootn-pair mesh
% Second subscript 2 indecates the second tootn-pair mesh

%============kh of healthy gear===========
kh=pi*E*L/(4*(1-v^2))*ones(length(kb11),1);
% kh=pi*E*L/(4*(1-v^2))*ones(roundn(2*pi/N1,-3)*1000,1);
%===================================
%
%=====Duration of double tooth-pair meshing=========
for theta1d=0.001:0.001:(roundn(thetaD,-3))
    %theta1d: angular displacement of the pinion
    % theta1s=thetaD:0.001:2*pi/N1;
    a11p=a110p+theta1d;                %for first pair of meshing teeth
    a11g=a110g-N1/N2*theta1d;    %for first pair of meshing teeth
    a12p=a11p+2*pi/N1;                 %for second pair of meshing teeth
    a12g=a11g-2*pi/N2;                 %for second pair of meshing teeth
    ir=roundn(theta1d,-3)*1000;
    %=======================================
    % Calculating ka, kb,ks of the First Tooth-Pair mesh of
    % the pinion and the gear during double tooth-pair meshing
    kb11(ir)=fkb(a11p,m,N1,E,L,a0);
    ka11(ir)=fka(a11p,m,N1,E,L,a0);
    ks11(ir)=fks(a11p,m,N1,E,L,v,a0);
    % vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    kb21(ir)=fkb(a11g,m,N2,E,L,a0);
    ka21(ir)=fka(a11g,m,N2,E,L,a0);
    ks21(ir)=fks(a11g,m,N2,E,L,v,a0);
    % Calculating ka, kb,ks of the Second Tooth-Pair mesh of
    % the pinion and the gear during Double tooth-pair meshing
    kb12(ir)=fkb(a12p,m,N1,E,L,a0);
    ka12(ir)=fka(a12p,m,N1,E,L,a0);
    ks12(ir)=fks(a12p,m,N1,E,L,v,a0);
    %=======================================
    kb22(ir)=fkb(a12g,m,N2,E,L,a0);
    ka22(ir)=fka(a12g,m,N2,E,L,a0);
    ks22(ir)=fks(a12g,m,N2,E,L,v,a0);
end
%=======================================
for theta1s=(thetaD+0.001):0.001:(roundn(2*pi/N1,-3)-0.001)
    % for theta1s=(thetaD+0.001):0.001:roundn(2*pi/N1,-3)
    a11p=a110p+theta1s;
    a11g=a110g-N1/N2*theta1s;
    a12p=a11p+2*pi/N1;
    a12g=a11g-2*pi/N2;
    jc=roundn(theta1s,-3)*1000;
    % Calculating ka, kb,ks of the First tooth-Pair mesh of
    % the pinion and the gear during Single tooth-pair meshing
    kb11(jc)=fkb(a11p,m,N1,E,L,a0);
    ka11(jc)=fka(a11p,m,N1,E,L,a0);
    ks11(jc)=fks(a11p,m,N1,E,L,v,a0);
    %=======================================
    kb21(jc)=fkb(a11g,m,N2,E,L,a0);
    ka21(jc)=fka(a11g,m,N2,E,L,a0);
    ks21(jc)=fks(a11g,m,N2,E,L,v,a0);
end
%===========Total effective meshing stiffness=========  总有效啮合刚度
kt1=1./(2./kh+1./kb11+1./ks11+1./ka11+1./kb21+1./ks21+1./ka21);
k2=1./(2./kh(1:size(kb12))+1./kb12+1./ks12+1./ka12+1./kb22+1./ks22+1./ka22);
kt2=zeros(length(kb11),1);
kt2(1:length(kb12))=k2;
kt=kt1+kt2;
plot(0.001*180/pi:0.001*180/pi:2*pi/N1*180/pi,kt,'k');
hold on
%=============connection line===================
angularD1=(0.001+pi*2/N1)*180/pi:0.001*180/pi:(roundn(thetaD,-3)+pi*2/N1)*180/pi;
kt01=kt(1:length(ka12));
plot(angularD1,kt01,'k');
angularS=(thetaD+0.001)*180/pi:0.001*180/pi:roundn(2*pi/N1,-3)*180/pi;
aa=[angularD1(1),kt01(1)];
bb=[angularS(end),kt(end)];
plot([aa(1),bb(1)],[aa(2),bb(2)],'k');
% %% 外啮合时变啮合刚度
% N1=17; N2=23; L=0.016; a0=20*pi/180; nv=75*pi/180;
% % nv denotes angle of start point of crack
% E=2.068e11; v=0.3; rintp=20; rintg=30;
% %===================================
% inva0=tan(a0)-a0;          a2p=pi/(2*N1)+inva0;
% inva0=tan(a0)-a0;          a2g=pi/(2*N2)+inva0;
% %===================================
% LOA=((N2+2)^2+(N1+N2)^2-2*(N2+2)*(N1+N2)*...
%     cos(acos(N2*cos(a0)/(N2+2))-a0))^0.5;
% a110p=tan(acos(N1*cos(a0)/LOA))-a2p;
% a110g=tan(acos(N2*cos(a0)/(N2+2)))-a2g;
% %====================================
% thetaD=tan(acos(N1*cos(a0)/(N1+2)))-2*pi/N1-...
%     tan(acos(N1*cos(a0)/LOA));
% %=======Initialization of matrix================
% kb11=zeros(roundn(2*pi/N1,-3)*1000-1,1);
% kb12=zeros(roundn(thetaD,-3)*1000,1);
% ka11=kb11; ka21=kb11;            ka12=kb12; ka22=kb12;
% ks11=kb11; ks21=kb11;            ks12=kb12; ks22=kb12;
% kb21=kb11;                               kb22=kb12;
% kf11=kb11; kf21=kb11;            kf12=kb12; kf22=kb12;
% % First subscript 1 or 2  indecates the pinion or the gear
% % Second subscript 1 indecates the first tootn-pair mesh
% % Second subscript 2 indecates the second tootn-pair mesh
% 
% %============kh of healthy gear===========
% kh=pi*E*L/(4*(1-v^2))*ones(length(kb11),1);
% % kh=pi*E*L/(4*(1-v^2))*ones(roundn(2*pi/N1,-3)*1000,1);
% %===================================
% %
% %=====Duration of double tooth-pair meshing=========
% for theta1d=0.001:0.001:(roundn(thetaD,-3))
%     %theta1d: angular displacement of the pinion
%     % theta1s=thetaD:0.001:2*pi/N1;
%     a11p=a110p+theta1d;                %for first pair of meshing teeth
%     a11g=a110g-N1/N2*theta1d;    %for first pair of meshing teeth
%     a12p=a11p+2*pi/N1;                 %for second pair of meshing teeth
%     a12g=a11g-2*pi/N2;                 %for second pair of meshing teeth
%     ir=roundn(theta1d,-3)*1000;
%     %=======================================
%     % Calculating ka, kb,ks of the First Tooth-Pair mesh of
%     % the pinion and the gear during double tooth-pair meshing
%     kb11(ir)=fkb(a11p,m,N1,E,L,a0);
%     ka11(ir)=fka(a11p,m,N1,E,L,a0);
%     ks11(ir)=fks(a11p,m,N1,E,L,v,a0);
%     % vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%     kb21(ir)=fkb(a11g,m,N2,E,L,a0);
%     ka21(ir)=fka(a11g,m,N2,E,L,a0);
%     ks21(ir)=fks(a11g,m,N2,E,L,v,a0);
%     % Calculating ka, kb,ks of the Second Tooth-Pair mesh of
%     % the pinion and the gear during Double tooth-pair meshing
%     kb12(ir)=fkb(a12p,m,N1,E,L,a0);
%     ka12(ir)=fka(a12p,m,N1,E,L,a0);
%     ks12(ir)=fks(a12p,m,N1,E,L,v,a0);
%     %=======================================
%     kb22(ir)=fkb(a12g,m,N2,E,L,a0);
%     ka22(ir)=fka(a12g,m,N2,E,L,a0);
%     ks22(ir)=fks(a12g,m,N2,E,L,v,a0);
% end
% %=======================================
% for theta1s=(thetaD+0.001):0.001:(roundn(2*pi/N1,-3)-0.001)
%     % for theta1s=(thetaD+0.001):0.001:roundn(2*pi/N1,-3)
%     a11p=a110p+theta1s;
%     a11g=a110g-N1/N2*theta1s;
%     a12p=a11p+2*pi/N1;
%     a12g=a11g-2*pi/N2;
%     jc=roundn(theta1s,-3)*1000;
%     % Calculating ka, kb,ks of the First tooth-Pair mesh of
%     % the pinion and the gear during Single tooth-pair meshing
%     kb11(jc)=fkb(a11p,m,N1,E,L,a0);
%     ka11(jc)=fka(a11p,m,N1,E,L,a0);
%     ks11(jc)=fks(a11p,m,N1,E,L,v,a0);
%     %=======================================
%     kb21(jc)=fkb(a11g,m,N2,E,L,a0);
%     ka21(jc)=fka(a11g,m,N2,E,L,a0);
%     ks21(jc)=fks(a11g,m,N2,E,L,v,a0);
% end
% %===========Total effective meshing stiffness=========  总有效啮合刚度
% kt1=1./(2./kh+1./kb11+1./ks11+1./ka11+1./kb21+1./ks21+1./ka21);
% k2=1./(2./kh(1:size(kb12))+1./kb12+1./ks12+1./ka12+1./kb22+1./ks22+1./ka22);
% kt2=zeros(length(kb11),1);
% kt2(1:length(kb12))=k2;
% kt=kt1+kt2;
% figure;
% plot(0.001*180/pi:0.001*180/pi:2*pi/N1*180/pi,kt,'k');
% hold on
% 
% 
% 
% %=============connection line===================
% angularD1=(0.001+pi*2/N1)*180/pi:0.001*180/pi:(roundn(thetaD,-3)+pi*2/N1)*180/pi;
% kt01=kt(1:length(ka12));
% plot(angularD1,kt01,'k');
% angularS=(thetaD+0.001)*180/pi:0.001*180/pi:roundn(2*pi/N1,-3)*180/pi;
% aa=[angularD1(1),kt01(1)];
% bb=[angularS(end),kt(end)];
% plot([aa(1),bb(1)],[aa(2),bb(2)],'k');
%% 刚度矩阵  psin为行星轮位置角 psin=2pi*（N-1)/N
kp1=1e8;
kp2=kp1;
kp3=kp2;
c33=[1e8 1e8 0]
psi1=0;
psi2=pi/3;
psi3=2*pi/3;
Kc11=kp1.*[1 0 -sin(psi1);0 1 cos(psi1);0 0 1];
Kc22=kp1.*[-cos(psi2) sin(psi2) 0;-sin(psi2) -cos(psi2) 0;0 -1 0];
Kc33=diag(c33);
phirp1=pi/9;
phirp2=140*pi/180;
phirp3=260*pi/180;
phisp1=-20*pi/180;
phisp2=100*pi/180;
phisp3=220*pi/180;
alphar=pi/9;
alphas=alphar;
alphap=alphar;
Ktt=2.296167998727640e+09; %内啮合平均啮合刚度
Ktt1=2.191207610444410e+09; %外啮合平均啮合刚度
Krp1=Ktt*[sin(phirp1)^2 -cos(phirp1)*sin(phirp1) -sin(phirp1);
    0 cos(phirp1)^2 cos(phirp1);
    0 0 1];
Krp2=Ktt*[-sin(phirp2)*sin(alphar) sin(phirp2)*cos(alphar) sin(phirp2);
    cos(phirp2)*sin(alphar) -cos(phirp2)*cos(alphar) -cos(phirp2);
    sin(alphar) -cos(alphar) -1];
Krp3=Ktt*[sin(alphar)^2 -cos(alphar)*sin(alphar) -sin(alphar);
    0 cos(alphar)^2 cos(alphar);
    0 0 1];
Ksp1=Ktt1*[sin(phisp1)^2 -cos(phisp1)*sin(phisp1) -sin(phisp1);
    0 cos(phisp1)^2 cos(phisp1);
    0 0 1];
Ksp2=Ktt1*[sin(phisp2)*sin(alphas) sin(phisp2)*cos(alphas) -sin(phisp2);
    -cos(phisp2)*sin(alphas) -cos(phisp2)*cos(alphas) cos(phisp2);
    -sin(alphas) -cos(alphas) 1];
Ksp3=Ktt1*[sin(alphas)^2 cos(alphas)*sin(alphas) -sin(alphas);
    0 cos(alphas)^2 -cos(alphas);
    0 0 1];
Kpp1=Kc11+Krp1+Ksp1;
Kpp2=Kc22+Krp2+Ksp2;
Kpp3=Kc33+Krp3+Ksp3;
TotalKcc=Kc11+Kc22+Kc33;
TotalKrp=Krp1+Krp2+Krp3;
TotalKsp=Ksp1+Ksp2+Ksp3;
Km=[TotalKcc o o Kc11 Kc22 Kc33;
    o TotalKrp o Krp1 Krp2 Krp3;
    o o TotalKsp Ksp1 Ksp2 Ksp3;
    o o o Kpp1 o o;
    o o o o Kpp2 o;
    o o o o o Kpp3];
%%  啮合误差
erp1=3e-7;
erp2=3e-7;
erp3=3e-7;
esp1=3e-7;
esp2=3e-7;
esp3=3e-7;
%%  激励 转速 负载等
O=[0 0 0];
Tc=100000;
Ts=28750;
Fr=Ktt*erp1*[sin(phirp1),-cos(phirp1),1];
Fs=Ktt1*esp1*[sin(phisp1),-cos(phisp1),1];
F1=Ktt*erp1*[sin(alphar),-sin(alphar),-1]+Ktt1*esp1*[sin(alphas),-sin(alphar),1];
Ft=[O Fr Fs F1 F1 F1]';
T=[0 0 -Tc/rc 0 0 0 0 0 Ts/rs 0 0 0 0 0 0 0 0 0]'
%%  固有频率
K=Kb+Km;
% 已知Kb Km G 
A=inv(M)*K;
[V,P] = eig(A);       %特征值和特征向量
%计算固有频率并排序
la=diag(P); %提取特征值
ww=sqrt(la);    %提取固有频率
w=sort(ww); %固有频率按阶次排序的向量
%提取特征值并按照阶次排序
N=length(M);
ki=0; %变换序列矩阵，模态向量{X}按阶次从小到大排序，所在列数字为原列数
    for j=1:N
    for i=1:N
    if w(j)==ww(i); X(:,j)=V(:,i)/max(abs(V(:,i)));ki(j)=i;end
    end
    end
    clc
%% 求位移 加速度 速度

%已知Kb18*18 M18*18 Komoge18*18 T18*1 Ft18*1 G18*18 行星架角速度172.5°/s
%q=[xc,yc,uc,xr,yr,ur,xs,ys,us,x1,y1,u1,x2,y2,u2,x3,y3,u3]; 18*18矩阵
%M*q''+wc*G*q'+[Kb+Km-wc^2*Komoge]*q==T+Ft




%% 啮合力Meshing Force
% MF=kt*0.005
% t=linspace(0,5,369);
% plot(t,MF,'k','linewidth',0.5);
% xlabel('时间 t/s')  %x轴坐标描述
% xlim([0 3])
% ylim([1.3 1.5])
% set(gca,'XTick',[0:1:5]);
% ylabel('接触力 F/N')
% %legend('外啮合副接触力');   %右上角标注
KKKKK=[2767043602.65021;2812612196.76322;2708658241.19066;2751976160.35579;2796708891.72790;2795520073.50203;2737290313.20951;2781214473.54904;2779785189.03421;2722973160.54521;2766114796.86873;2810651375.20011;2809640364.88559;2751396370.64110;2795129492.90400;2793875226.75577;2737046319.62391;2780001820.47172;2778518927.39209;2723052348.59503;2765254901.71062;2808787087.15571;2807697817.20631;2750875893.48373;2793635652.02710;2792315260.32921;2736852530.14969;2778864322.65446;2822185938.43250;2723173090.48309;2764460294.61483;2807014961.88466;2805847292.24946;2750411338.08660;2792223293.25880;2790836001.59950;2736705764.99839;2777798173.36773;2820138669.69251;2723332398.43698;2763727414.36883;2805330716.71597;2804084388.16688;2749999367.16588;2790888417.85103;2789433366.78350;2736602891.04818;2776799633.30045;2818185738.93507;2723527325.38783;2763052758.67891;2803730153.30117;2802404800.46515;2749636694.90511;2789627100.95068;2788103357.11440;2736540820.13065;2775865027.84775;2816322740.92822;2723754963.54534;2762432881.87762;2802209154.07090;2800804317.74650;2749320085.02678;2788435488.58245;2786842055.23310;2736516507.41851;2774990744.55155;2814545358.51010;2813250408.95959;2761864392.76210;2800763678.88979;2799278817.69061;2749046348.97655;2787309794.80153;2785645621.78293;2736526949.91027;2774173230.68674;2812849358.89184;2811464312.97382;2761343952.55778;2799389761.89966;2797824263.26423;2748812344.21534;2786246299.00816;2824692508.37586;2736569185.00806;2773408990.98668;2811230590.16933;2809755271.95603;2760868273.00191;2798083508.54231;2796436699.14632;2748614972.61493;2785241343.41659;2822832436.19486;2736640289.18492;2772694585.50159;2809684978.03350;2808119133.58661]
%KKKKK是截取的一段时变啮合刚度
MF=KKKKK*0.000005; %0.000005是变形
t=linspace(0,0.3,100);
figure
plot(t,MF,'k','linewidth',0.5);



