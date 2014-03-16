% Exercício prático sobre teoria de Ekman.
%

% Limpando tudo
clear all; close all; clc;

% Lendo dados de vento.
load('AS_vento.mat');
[long,latg] = meshgrid(lon,lat);
lon2 = lon(2:end-1);
lat2 = lat(2:end-1);
long2 = long(2:end-1,2:end-1);
latg2 = latg(2:end-1,2:end-1);

%+++++++++++++++++++++++++++++++++++++++++++++
% Parte A - A camada de Ekman de superfície.
%           Transporte e bombeamento de Ekman.

%==================== Calculando o vetor tensão de cisalhamento do vento. ====================%
% Verão.
rho_ar = 1.226;                     % kg/m3
magver = sqrt(uver.^2 + vver.^2);   % m/s
Cd = 7.5e-4 + 6.7e-5*magver;        % Mellor (2004).
tauxver = rho_ar.*Cd.*magver.*uver; % N/m2
tauyver = rho_ar.*Cd.*magver.*vver; % N/m2

% Inverno.
rho_ar = 1.226;                     % kg/m3
maginv = sqrt(uinv.^2 + vinv.^2);   % m/s
Cd = 7.5e-4 + 6.7e-5*maginv;        % Mellor (2004).
tauxinv = rho_ar.*Cd.*maginv.*uinv; % N/m2
tauyinv = rho_ar.*Cd.*maginv.*vinv; % N/m2

%========================== Calculando o vetor transporte de Ekman. ==========================%
rho0 = 1025.;                       % kg/m3, densidade média da água do mar na superfície.
omega = 2*pi/86400.;                % rad/s
f0 = 2.*omega.*sind(latg);          % rad/s

% Verão.
Uever = +tauyver ./ (rho0.*f0);     % m2/s
Vever = -tauxver ./ (rho0.*f0);     % m2/s

% Inverno.
Ueinv = +tauyinv ./ (rho0.*f0);     % m2/s
Veinv = -tauxinv ./ (rho0.*f0);     % m2/s

%============================= Calculando o bombeamento de Ekman. ============================%
%=========== é a velocidade vertical forçada pelo vento na base da camada de Ekman. ==========%

% Matrizes de distância meridional (dy) e zonal (dx).
lat2m = 60.*1852.;         % [m]
d = 0.25;                  % A resolução é de 1/4 de grau.
dy = d.*lat2m;             % A distância de 1 grau de latitude é constante.
dx = d.*lat2m.*cosd(latg); % A distância de 1 grau de longitude depende da latitude.
dx = dx(2:end-1,2:end-1);
f0 = f0(2:end-1,2:end-1); % Cortando as matrizes para ficarem do mesmo tamanho que vx e uy.

% Verão.
dvdx = +(tauyver(2:end-1,3:end) - tauyver(2:end-1,1:end-2)) ./ (2.*dx);
dudy = -(tauxver(3:end,2:end-1) - tauxver(1:end-2,2:end-1)) ./ (2.*dy); % O menos é pq o vetor lat está ao contrário.

rot = dvdx - dudy;         % [N/m3] componente vertical do rotacional do vetor tensão de cisalhamento do vento.

Wever = rot ./ (rho0.*f0); % [m/s] Bombeamento de Ekman.

% Eliminando valores inconsistentes.
msk = isinf(Wever);
Wever(msk) = NaN;

% Inverno.
dvdx = +(tauyinv(2:end-1,3:end) - tauyinv(2:end-1,1:end-2)) ./ (2.*dx);
dudy = -(tauxinv(3:end,2:end-1) - tauxinv(1:end-2,2:end-1)) ./ (2.*dy); % O menos é pq o vetor lat está ao contrário.

rot = dvdx - dudy;         % [N/m3] componente vertical do rotacional do vetor tensão de cisalhamento do vento.

Weinv = rot ./ (rho0.*f0); % [m/s] Bombeamento de Ekman.

% Eliminando valores inconsistentes.
msk = isinf(Weinv);
Weinv(msk) = NaN;

%=============================================================================================%

xp = 0; % Longitude do perfil meridional que vai ser plotado.
Xp = xp.*ones(length(lat),1);

% 1) Campos climatológicos de verão (janeiro).
grey1 = [.95,.95,.95];
grey2 = [.6,.6,.6];

figure;
set(gcf,'color',grey1);
subplot(131)
m_proj('mercator','lon',[-60 20],'lat',[-60 -12],'on');
hold on
set(gcf,'color','w')
ds = 91; % Quantos vetores pular para cada vetor plotado.
scale = 2.0;
m_quiver(long(1:ds:end),latg(1:ds:end),tauxver(1:ds:end),tauyver(1:ds:end),scale,'color','r');
m_usercoast('dados/coast.mat','patch',grey1,'LineStyle','-');
m_plot(Xp,lat,'k--','linewidth',1.8); % Marca o perfil plotado ao lado.
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',10,'linewidth',0);
title('Tensão de cisalhamento - janeiro','fontsize',14,'fontweight','bold');
hold off

subplot(132)
hold on
plot(Vever(:,lon==xp),lat,'color',grey2,'linewidth',2.0);
% plot(tauxver(:,lon==xp),lat,'color','g','linewidth',2.0);
plot([0,0],[min(lat),max(lat)],'k--','linewidth',1.8);
pbaspect([1,1,1]);
xlabel('V_e [m^2 s^{-1}]','fontsize',12,'fontweight','bold');
ylabel('Latitude','fontsize',12,'fontweight','bold');
xlim([-2,2]);
grid on;
title('Componente meridional do transporte','fontsize',14,'fontweight','bold');
hold off

subplot(133)
hold on
plot(Wever(:,lon2==xp),lat2,'m-','linewidth',2.0);
plot([0,0],[min(lat),max(lat)],'k--','linewidth',1.8);
pbaspect([1,1,1]);
xlabel('W_e [m s^{-1}]','fontsize',12,'fontweight','bold');
ylabel('Latitude','fontsize',12,'fontweight','bold');
xlim([-0.5e-5,0.5e-5]);
grid on;
title('Bombeamento de Ekman','fontsize',14,'fontweight','bold');
hold off

clc
disp('Vento climatológico de verão.')

pause;

% 2) Campos climatológicos de inverno (julho).

figure;
set(gcf,'color',grey1);
subplot(131)
m_proj('mercator','lon',[-60 20],'lat',[-60 -12],'on');
hold on
set(gcf,'color','w')
ds = 91; % Quantos vetores pular para cada vetor plotado.
scale = 2.0;
m_quiver(long(1:ds:end),latg(1:ds:end),tauxinv(1:ds:end),tauyinv(1:ds:end),scale,'color','r');
m_usercoast('dados/coast.mat','patch',grey1,'LineStyle','-');
m_plot(Xp,lat,'k--','linewidth',1.8); % Marca o perfil plotado ao lado.
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',10,'linewidth',0);
title('Tensão de cisalhamento - julho','fontsize',14,'fontweight','bold');
hold off

subplot(132)
hold on
plot(Veinv(:,lon==xp),lat,'color',grey2,'linewidth',2.0);
plot([0,0],[min(lat),max(lat)],'k--','linewidth',1.8);
pbaspect([1,1,1]);
xlabel('V_e [m^2 s^{-1}]','fontsize',12,'fontweight','bold');
ylabel('Latitude','fontsize',12,'fontweight','bold');
xlim([-2,2]);
grid on;
title('Componente meridional do transporte','fontsize',14,'fontweight','bold');
hold off

subplot(133)
hold on
plot(Weinv(:,lon2==xp),lat2,'m-','linewidth',2.0);
plot([0,0],[min(lat),max(lat)],'k--','linewidth',1.8);
pbaspect([1,1,1]);
xlabel('W_e [m s^{-1}]','fontsize',12,'fontweight','bold');
ylabel('Latitude','fontsize',12,'fontweight','bold');
xlim([-0.5e-5,0.5e-5]);
grid on;
title('Bombeamento de Ekman','fontsize',14,'fontweight','bold');
hold off

clc
disp('Vento climatológico de inverno.')
disp('Exercício: Interprete e compare estas 2 figuras. Pense em convergência e divergência.')

pause;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Parte B - A camada de Ekman de superfície.
%           O efeito do Av variável na espiral de Ekman.

lat=-25;
f0=sw_f(lat);
af0=abs(f0);
s=f0/af0;

% Parâmetros.

Av=0.045;           % Valor característico de Av.
h=sqrt(2*Av/af0);   % espessura da camada de Ekman.
he=pi*h;            % espessura efetiva da camada de Ekman.
H=100;              % máxima profundidade do modelo.
dz=.5;              % incremento de profundidade.
z=-H:dz:0;          % Eixo vertical.
del=s*af0*dz*dz;

nmax=H/dz + 1;

% Aproximações teóricas para o perfil de Av(z)

% constante.
Avc=Av*ones(nmax,1);
% linear.
Avl = (40 + .39.*z).*1e-3;
% exponencial.
Ave = (2.5 + 40.*exp(z/10.5))*1e-3;

% Vento.
uwind=10; vwind=0;
wind=sqrt(uwind^2+vwind^2);
rhoa=1.22;
CD=1e-3*(0.8+0.065*wind);  % from Wu(JGR,1982)
taux=rho_ar*CD*wind*uwind;
tauy=rho_ar*CD*wind*vwind;
Tau=taux+i*tauy;

%+++ Resolvendo o problema numericamente com a técnica de mínimos quadrados.

% Construindo as matrizes tridiagonais A.
Avcdel=del./Avc;
Avldel=del./Avl;
Avedel=del./Ave;
cA= diag(2+i*Avcdel) + diag(-ones(nmax-1,1),1) + diag(-ones(nmax-1,1),-1); % constante.
lA= diag(2+i*Avldel) + diag(-ones(nmax-1,1),1) + diag(-ones(nmax-1,1),-1); % linear.
eA= diag(2+i*Avedel) + diag(-ones(nmax-1,1),1) + diag(-ones(nmax-1,1),-1); % exponencial.

% Impondo as condições de contorno.
% No fundo
cA(1,2)=0; cA(1,1)=1;
lA(1,2)=0; lA(1,1)=1;
eA(1,2)=0; eA(1,1)=1;
% Na superfície.
cA(nmax,nmax)= 1;
lA(nmax,nmax)= 1;
eA(nmax,nmax)= 1;
% Criando os vetores resposta D.
cD=zeros(nmax,1);
lD=zeros(nmax,1);
eD=zeros(nmax,1);
% Consertando a condição de contorno superior nos vetores D.

cD(nmax)= dz/(rho0*Avl(nmax))*Tau;
lD(nmax)= dz/(rho0*Avl(nmax))*Tau;
eD(nmax)= dz/(rho0*Ave(nmax))*Tau;

% Solução do problema: A*V=D;

Vc=cA\cD;
Vl=lA\lD;
Ve=eA\eD;

u_c=real(Vc); v_c=imag(Vc);
u_l=real(Vl); v_l=imag(Vl);
u_e=real(Ve); v_e=imag(Ve);

%%% Plotando.

% Perfis de Av.
asp = [1,1,1];
figure;
ax1 = subplot(1,3,1);
plot(Avc,z,'b',Avl,z,'r',Ave,z,'k','linewidth',2);
hold on
grid on;
hold off;
xlabel('Av [m^2/s]','fontsize',12,'fontweight','bold')
ylabel('Profundidade [m]','fontsize',12,'fontweight','bold')
title('Perfis do coeficiente de viscosidade turbulenta','fontsize',12,'fontweight','bold')
pbaspect([asp]);

% Vetores.
ax2 = subplot(1,3,2);
compass(Ve(1:10:end),'k')
hold on
compass(Vl(1:10:end),'r')
compass(Vc(1:10:end),'b')
hold off
title('Vetores de velocidade','fontsize',12,'fontweight','bold')
pbaspect([asp]);

% Hodógrafos dos vetores.
ax3 = subplot(1,3,3);
plot(u_c,v_c,'b',u_l,v_l,'r',u_e,v_e,'k','linewidth',1.7);
grid on;
xlabel('u [m/s]','fontsize',12,'fontweight','bold')
ylabel('v [m/s]','fontsize',12,'fontweight','bold')
title('Hodógrafos dos vetores','fontsize',12,'fontweight','bold')
pbaspect([asp]);

axes(ax3);
legend('A_v constante','A_v linear','A_v exponencial','location','southeast')

clc
disp('Espiral de Ekman superficial com Av = Av(z).')
disp('Exercício: Calcule o ângulo entre o vento e a corrente na superfície para cada caso de Av(z).')

pause;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Parte C - A camada de Ekman de fundo.
%           Visualização da espiral de Ekman de fundo.

% Parâmetros.
lat=-25;
f0=sw_f(lat);
af0=abs(f0);
s=f0/af0;
Av=0.045;           % Valor característico de Av.
h=sqrt(2*Av/af0);   % espessura da camada de Ekman.
he=pi*h;            % espessura efetiva da camada de Ekman.
H=300;              % Altura máxima sobre o fundo.
dz=.5;              % incremento de profundidade.
z=0:dz:H;           % Eixo vertical.
del=s*af0*dz*dz;

% Componentes do vetor corrente geostrófica.
ug = 0.9; % [m/s]
vg = 0;   % [m/s]

% Soluções analíticas, para Av constante.
% e.g., Vallis (2006).
u = ug - exp(-z/h) .* ( ug.*cos(z/h) + vg.*sin(z/h) );
v = vg + exp(-z/h) .* ( ug.*sin(z/h) - vg.*cos(z/h) );

V = u + i.*v;

%%% Plotando.

% Perfis de Av.
asp = [1,1,1];
figure;
plot(u,z,'b',v,z,'r','linewidth',2);
hold on
xl = xlim;
plot(xl,[h h],'k--','linewidth',2);
plot(xl,[he he],'color',[.5 .5 .5],'linestyle','--','linewidth',2);
grid on;
xt = 0.1;yt = 270;
text(xt,yt,['ug = ',num2str(ug),' m/s'],'fontsize',15,'fontweight','bold')
text(xt,yt-20,['vg = ',num2str(vg),' m/s'],'fontsize',15,'fontweight','bold')
legend('u','v','location','southeast')
hold off;
xlabel('u,v [m/s]','fontsize',15,'fontweight','bold')
ylabel('Altura acima do fundo [m]','fontsize',15,'fontweight','bold')
title('Perfis de velocidade total na camada de Ekman de fundo.','fontsize',12,'fontweight','bold')
pbaspect([asp]);

% Vetores.
figure;
compass(V(1:30:end),'m')
title('Vetor velocidade','fontsize',12,'fontweight','bold')
pbaspect([asp]);

% Hodógrafos dos vetores.
figure;
plot(u,v,'m-','linewidth',2.0);
grid on;
xlabel('u [m/s]','fontsize',12,'fontweight','bold')
ylabel('v [m/s]','fontsize',12,'fontweight','bold')
title('Hodógrafo do vetor velocidade','fontsize',12,'fontweight','bold')
pbaspect([asp]);

clc;
disp('O quê as linhas tracejadas marcam?')
pause;
