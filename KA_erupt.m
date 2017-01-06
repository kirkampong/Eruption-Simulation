% Steady state conduit flow model
% Kirk Ampong
% June 2016

% user input parameters
K.opt1 = input('Enter "constant" or "dynamic" for preferred viscosity model:');
K.opt2 = input('Enter "standard" or "improved" to choose crystal viscosity effect model: '); 
K.opt3 = input('Enter "shear" or "gas_fraction" to choose fragmentation model: ');

% properties:

H = 8e3; % conduit length                    
rhor = 2500; % rock density, kg/m^3
K.r = 30; % conduit radius, m
K.X0 = 0.04; % total mass concentration of volatiles
%K.n0 = 0.04; % total mass fraction of volatiles

% solubility of water in basalt
%A.s = 6.3e-8; % solubility constant (for p in Pa)
%A.m = 0.7; % solubility law exponent

% solubility of water in rhyolite
K.s = 4e-6; % solubility constant (for p in Pa)
K.m = 0.5; % solubility law exponent

K.R = 461; % specific gas constant of volatiles, J/kgK
K.T = 1000; % temperature, K

% liquid
K.rhol0 = 2600; % reference liquid density, kg/m^3
K.pl0 = 10e6; % reference pressure (for liquid density), Pa
K.Kl = 10e9; % liquid bulk modulus, Pa

% crystals
K.PHI = 0;  % crystal mass fraction
K.rhoc0 = 2600; % reference crystal density, kg/m^3
K.Kl_c = 10e9; % crystal bulk modulus, Pa

K.f0 = 0.01; % Darcy-Weisbach friction factor
K.phi0 = 0.75; % critical gas volume fraction for fragmentation

K.g = 9.8; % acceleration due to gravity, m/s^2

% pressure at chamber

Deltap = 1e6; % excess pressure, Pa
p0 = rhor*K.g*H+Deltap; % lithostatic + excess pressure


K.n0 = K.X0.*(1-K.PHI)./(1+K.X0); %total mass fraction of volatiles
%K.X0 = K.n0/(1-K.n0-K.PHI); % total mass concentration of volatiles

% solve using shooting method by varying initial velocity
% u0 to obtain choked flow at vent
% shooting method implemented using fzero
warning off MATLAB:ode15s:IntegrationTolNotMet
opt = odeset('NormControl','on','RelTol',1e-12,'AbsTol',1e-12);
u0 = fzero(@(u0) KA_topBC(u0,p0,K,H,opt), 1, optimset('Display','iter'));

% evaluate solution
u0 = 0.999*u0; % slightly below velocity for choked flow at vent
sol = ode15s(@KA_eruptODE,[-H 0],[p0; u0],opt,K); % solve ODE system
y = sol.x;
[q,dqdy] = deval(sol,y);
p = q(1,:); u = q(2,:);
dpdy = dqdy(1,:); dudy = dqdy(2,:);
[rho,phi,c,beta,rhoc,rhol,rhoe,rhod,Xd,Xe,ne,nd] = KA_eos(p,K); % density and other fields
dudy = -beta.*u.*dpdy;
tau = KA_wallshear(rho,phi,rhoc,rhol,rhoe,rhod,Xd,Xe,u,K,c,beta); % wall shear stress, Pa
Q = rho(1)*u(1)*pi*K.r^2; % mass eruption rate, kg/s


%-------------------------------------------------------------------------

% alternative units for plotting (Y in km, u in m/s, P in MPa)
% Y = y*1e-3; P = p*1e-6;
% plot density vs depth
% figure, plot(rho,Y), xlabel('density, (kg/m^3)'),ylabel('depth, (km)')
% plot velocity vs depth
% figure, plot(u,Y), xlabel('velocity, (m/s)'),ylabel('depth, (km)')

% add lines showing exsolution (black) and fragmentation (red) depths
% Yex = Y(find(Xe>0,1)); % exsolution depth
% Yfr = Y(find(phi>K.phi0,1)); % fragmentation depth for gas fraction model
% xl=xlim;
% hold on
% plot([min(xl) max(xl)],[Yex Yex],'k--')
% plot([min(xl) max(xl)],[Yfr Yfr],'r--')
% hold off

