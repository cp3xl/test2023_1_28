

clc
clear
%% Compressor parameters
x=[0.000459/5.81,20.5,0.00969,0.01025,28,80,293,0.874 70000000];
A1=x(1);   %0.00015;
u=x(2);   %气体动力粘度
r2=x(3);   %0.015 ;
r1=x(4);   %0.01 
B1=x(5)*pi/180;   % 41
B2=x(6)*pi/180;   % 80
T01=x(7);
c=x(8)*0.92; % 
Kf=x(9);
Cp=   1006  ;       %恒压下的具体温度
K=    1.4   ;     %K=Cp/Cv(Cv是恒体积下的具体温度)
q1=   1.29  ;      %入口密度
D1=2*r1;
Li=0.044;%假设的平均叶轮通道长度
Ld=0.040;%假设的平均扩散器通道长度
Di=0.018;%假设的平均叶轮液压通道直径
Dd=0.017;%假设的平均扩散器液压通道直径

%% High altitude 2000m atmosphere

h=0;
T0=293.0;
P0=101325;    % Pa
L=-0.0065;    % K/m
g0=9.80665;
Ma=0.0289644;
R=8.31432;

if (h<11000)
Ph=P0*(1+(L/T0)*(h-0))^(-g0*Ma/R/L);
Th=T0+h*L;
else if (11000<=h<=20000)
Ps=P0*(1+(L/T0)*(11000-0))^(-g0*Ma/R/L);
Th=T0+11000*L;
Ph=Ps*exp(-g0*Ma*(h-11000)/R/Th);
    else
        Th=0;
        Ph=0;
    end
end

q1=Ph*Ma/R/Th*(1-0.0061*(1-18.0/28.97));

a=(cot(B1))/(q1*A1*r1);
D=(cot(B2))/(q1*A1*r1);
G1=(c*(r2)^2-0.5*(r1)^2)/(Cp*Th) ;
G2=-(0.5*(r1)^2*a^2)/(Cp*Th) ;
G3=((r1)^2*a-(r2)^2*c*D)/(Cp*Th) ;

%% Compressor map
n(1)=100000;
d=100000*ones(1,26);
w(1)=(pi/30.0)*n(1);
m=0:0.001:0.025;

for i=1:length(m)
    Re(1)=D1*D1*pi*q1/(60*g0*u*10^(-6))*n(1);
    f(1)=0.3164*(Re(1))^(-0.25);
    kf(1)=2*f(1)*Li/(Di*q1*q1*A1*A1*sin(B1)*sin(B1))+2*f(1)*Ld/(Dd*q1*q1*A1*A1*sin(B2)*sin(B2));
    G4(1)=kf(1);
P(i)=Ph*(1+G1*w(1)*w(1)+G2*m(i)*m(i)+G3*w(1)*m(i)+G4(1)*m(i)*m(i)/(-Cp*Th))^(K/(K-1));


end


figure(1);
plot3(m,P,d);

figure(2);
[m,P]=meshgrid(m,P);
mesh(m,d,P);

%axis([0 0.025 Ph 1.1*max(P5)])
%axis([0 0.025 Ph 40000])
%xlabel('speed (kg/s)','FontSize',14);
%ylabel('Air flow (kg/s)','FontSize',14);
%zlabel('Air pressure (Pa)','FontSize',14);
%set(gca,'FontSize',14)
grid on;


%legend('150,000 rpm','180,000 rpm','210,000 rpm','240,000 rpm','280,000 rpm');
%gtext('h=2000 m','FontSize',14);



