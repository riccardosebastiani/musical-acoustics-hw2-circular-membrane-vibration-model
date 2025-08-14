%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            HW 2                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
%% DATA

r_m = 0.15;   % [m]
T = 10;     % [N/m]
sigma = 0.07;   % [kg/m^2]  unit surface weight (area mass density)

r_struck = 0.075;    %[m]
phi_struck = 15*(2*pi)/360;   %[rad]

r_sensor = 0.075;    %[m]
phi_sensor = 195*(2*pi)/360;   %[rad]

Q = 25;

t = 0:0.001:0.1;
force = 0.1*exp(-(t-0.03).^2./0.01^2);

%% a) Propagation speed in the membrane

v = sqrt(T/sigma);  % for transverse waves in a clamped edge membrane

%% b) Frequency of the first eighteen modes

rSpatialNPoints = 100;
phiSpatialNPoints = 100;
A = 1;

dc = 0.1;
c = 0:dc:20;
J = zeros(6,length(c));
for i = 0:5
    J(i+1,:) = besselj(i,c);
end


% % figure()
% % plot(z,J)
% % grid on
% % legend('J_0','J_1','J_2','J_3','J_4','J_5','Location','Best')
% % title('Bessel Functions of the First Kind for $\nu \in [0, 5]$','interpreter','latex')
% % xlabel('z','interpreter','latex')
% % ylabel('$J_\nu(z)$','interpreter','latex')


%from slides
J_zeros = [2.4048 3.8317 5.1356 6.3802 7.5883 8.7715 9.936 11.086 12.225;
           5.5201 7.0156 8.4172 9.7610 11.0647 12.3386 50 50 50;
           8.6357 10.1735 11.6198 13.0152 14.3725 15.7002 50 50 50;
           11.7915 13.3237 14.7960 16.2235 17.6160 18.9801 50 50 50;
           14.9309 16.4706 17.9598 19.4094 20.8269 22.2178 50 50 50];

f_mn = v/(2*pi*r_m)*J_zeros;
f_linear = reshape(f_mn,1,[]);  %f_mn into linear array
f_linear = mink(f_linear,18);   %first 18 smallest elements of f_linear
omega_linear = 2*pi*f_linear;


% % theta = linspace(0,2*pi,300) ; 
% % r = 0:0.1:1 ; 
% % [T,R] = meshgrid(theta,r) ; 
% % X = R.*cos(T) ;
% % Y = R.*sin(T) ; 
% % Z = sqrt(X.^2+Y.^2) ;
% % %surf(X,Y,Z)
% % pcolor(X,Y,Z)
% % shading interp
% % axis equal

radius = linspace(0,r_m,rSpatialNPoints);
phi = linspace(0,2*pi,phiSpatialNPoints);

z = zeros(length(phi),length(radius));
z_modes = zeros(length(phi),length(radius),18);

zero = 3.8317;
zero_vec = reshape(J_zeros,1,[]);
zero_vec_min = mink(zero_vec,18);


x = zero*radius/r_m;
Jm(1,:) = besselj(2,x);

x_vec = zero_vec_min'*radius/r_m;

Jm_vec = zeros(length(zero_vec_min),length(x));

Jm_vec(1,:) = besselj(1,x_vec(1,:));    %m=0, n=1
Jm_vec(2,:) = besselj(2,x_vec(2,:));    %m=1, n=1
Jm_vec(3,:) = besselj(3,x_vec(3,:));    %m=2, n=1
Jm_vec(4,:) = besselj(1,x_vec(4,:));    %m=0, n=2
Jm_vec(5,:) = besselj(4,x_vec(5,:));    %m=3, n=1
Jm_vec(6,:) = besselj(2,x_vec(6,:));    %m=1, n=2
Jm_vec(7,:) = besselj(5,x_vec(7,:));    %m=4, n=1
Jm_vec(8,:) = besselj(3,x_vec(8,:));    %m=2, n=2
Jm_vec(9,:) = besselj(1,x_vec(9,:));    %m=0, n=3
Jm_vec(10,:) = besselj(6,x_vec(10,:));  %m=5, n=1
Jm_vec(11,:) = besselj(4,x_vec(11,:));  %m=3, n=2
Jm_vec(12,:) = besselj(7,x_vec(12,:));  %m=6, n=1
Jm_vec(13,:) = besselj(2,x_vec(13,:));  %m=1, n=3
Jm_vec(14,:) = besselj(5,x_vec(14,:));  %m=4, n=2
Jm_vec(15,:) = besselj(8,x_vec(15,:));  %m=7, n=1
Jm_vec(16,:) = besselj(3,x_vec(16,:));  %m=2, n=3
Jm_vec(17,:) = besselj(1,x_vec(17,:));  %m=0, n=4
Jm_vec(18,:) = besselj(9,x_vec(18,:));  %m=8, n=1


for i = 1:length(x)
    
    z(:,i) = A * real(exp(1i*phi)) * Jm(1,i);
    
end

z_tr = z.';


for j=1:length(x)
        z_modes(:,j,1) =  A * real(exp(1i*0*phi)) .* Jm_vec(1,j);
        z_modes(:,j,2) =  A * real(exp(1i*phi)) .* Jm_vec(2,j);
        z_modes(:,j,3) =  A * real(exp(1i*2*phi)) .* Jm_vec(3,j);
        z_modes(:,j,4) =  A * real(exp(1i*0*phi)) .* Jm_vec(4,j);
        z_modes(:,j,5) =  A * real(exp(1i*3*phi)) .* Jm_vec(5,j);
        z_modes(:,j,6) =  A * real(exp(1i*phi)) .* Jm_vec(6,j);
        z_modes(:,j,7) =  A * real(exp(1i*4*phi)) .* Jm_vec(7,j);
        z_modes(:,j,8) =  A * real(exp(1i*2*phi)) .* Jm_vec(8,j);
        z_modes(:,j,9) =  A * real(exp(1i*0*phi)) .* Jm_vec(9,j);
        z_modes(:,j,10) =  A * real(exp(1i*5*phi)) .* Jm_vec(10,j);
        z_modes(:,j,11) =  A * real(exp(1i*3*phi)) .* Jm_vec(11,j);
        z_modes(:,j,12) =  A * real(exp(1i*6*phi)) .* Jm_vec(12,j);
        z_modes(:,j,13) =  A * real(exp(1i*phi)) .* Jm_vec(13,j);
        z_modes(:,j,14) =  A * real(exp(1i*4*phi)) .* Jm_vec(14,j);
        z_modes(:,j,15) =  A * real(exp(1i*7*phi)) .* Jm_vec(15,j);
        z_modes(:,j,16) =  A * real(exp(1i*2*phi)) .* Jm_vec(16,j);
        z_modes(:,j,17) =  A * real(exp(1i*0*phi)) .* Jm_vec(17,j);
        z_modes(:,j,18) =  A * real(exp(1i*8*phi)) .* Jm_vec(18,j);
end

[T,R] = meshgrid(phi,radius) ; 
X = R.*cos(T) ;
Y = R.*sin(T) ;  

% % figure,
% % for i =1:size(z_tr,2)
% %     surf(X,Y,z_tr)
% %     xlabel('$r \cos \phi$','Interpreter','latex','Fontsize',15)
% %     ylabel('$r \sin\phi$','Interpreter','latex', 'Fontsize',15)
% %     zlabel('$z_{m n}$','Interpreter','latex','Fontsize',15)
% %     zlim([-A A])
% %     colorbar
% %     title('Mode (1,2)','Fontsize',20)
% %     %xlim([0 x(end)])
% %     %ylim([0 theta(end)])
% % end


%% c) displacement time signal read by sensor

figure()
plot(t,force,'LineWidth',1.2)

% % a=besselj(0, z) - besselj(1, z)./z;
% % dJn = zeros(5,length(z));
% % for i = 1:4
% %     dJn(i,:) = J(i,:) - J(i+1,:)./z;
% % end
% % 
% % figure()
% % plot(z,dJn)
% % grid on
% % legend('dJ_0','dJ_1','dJ_2','dJ_3','dJ_4','Location','Best')
% % title('Derivative of Bessel Functions of the First Kind for $\nu \in [0, 4]$','interpreter','latex')
% % xlabel('z','interpreter','latex')
% % ylabel('$dJ_\nu(z)$','interpreter','latex')


% % dJ = zeros(6,length(z));
% % for i = 1:6
% %     dJ(i,:) = gradient(J(i,:),z);  % or: dJ = gradient(J)./gradient(z);
% % end
% % 
% % 
% % figure()
% % plot(z,dJ)
% % grid on
% % legend('dJ_0','dJ_1','dJ_2','dJ_3','dJ_4','dJ_5','Location','Best')
% % title('Derivative of Bessel Functions of the First Kind for $\nu \in [0, 5]$','interpreter','latex')
% % xlabel('z','interpreter','latex')
% % ylabel('$dJ_\nu(z)$','interpreter','latex')


mass = zeros(1,18);
stiffness = zeros(1,18);
damping = zeros(1,18);


for i=1:18
    fun = @(radius) ((sigma*2*pi*radius) .* abs(z_modes(:,:,i)).^2);
    q = trapz(radius,fun(radius),2);
    mass(i) = q(1);
    stiffness(i) = omega_linear(i)^2 * mass(i); 
    %stiffness(i) = omega_linear(i)^2;
    damping(i) = sqrt(mass(i)*stiffness(i))/Q;
    %damping(i) = stiffness(i)/Q;
end



%receptance

x_vec_struck = zero_vec_min*r_struck/r_m;
x_vec_sensor = x_vec_struck;

m = [0 1 2 0 3 1 4 2 0 5 3 6 1 4 7 2 0 8];
Jm_vec_struck = zeros(length(m),1);
num = zeros(1,length(m));

for i = 1:length(m)
    Jm_vec_struck(i) = besselj(m(i),x_vec_struck(i));
end


for i = 1:length(m)
    num(i) = A * real(exp(1i*m(i)*(phi_struck+phi_sensor))) * Jm_vec_struck(i)^2;
end


w = linspace(0,omega_linear(end),3000);
H = zeros(length(w),length(m));
H_sum = zeros(1,length(w));

for i=1:length(m)
    for ii=1:length(w)
        H(ii,i) = (num(i))/(-w(ii)^2 * mass(i) + 1i*w(ii)*damping(i) + stiffness(i));
        H_sum(ii) = H_sum(ii) + H(ii,i);
    end
end

h_sum = ifft(H_sum);



%c = conv(h_sum,force);

%time = linspace(0,0.0001,length(c));

time_h = linspace(0,0.1,length(H_sum));

force_h = 0.1*exp(-(time_h-0.03).^2./0.01^2);

c = conv(h_sum,force_h);


force_w = fft(force_h);
z_w = force_w.*H_sum;
z_t = ifft(z_w);

t_plot = linspace(0,0.1,length(c));

figure()
% plot(time,real(c))
% hold on
% plot(time,imag(c))
% hold on
plot(time_h,force_h)
title('force')
xlabel('time [s]')
ylabel('f(t) [N]')
    
figure()
% plot(time,real(c))
% hold on
% plot(time,imag(c))
% hold on
plot(t_plot, real(c))
ylabel('z(t) [m]')
xlabel('time [s]')
title('displacement')


%% DATA ex2

h = 0.001;  %[m] thickness
E_p = 69e9; %[Pa];
rho_p = 2700; %[kg/m^3];
ni = 0.334;

V_p = h*pi*r_m^2;
mass_p = rho_p*V_p;

%% d) propagation speed quasi-longitudinal and longitudinal waves

c_ql = sqrt(E_p/(rho_p*(1-ni^2)));

c_l = sqrt((E_p*(1-ni))/(rho_p*(1+ni)*(1-2*ni)));


%% e) propagation speed bending waves

f = 0:0.1:2*f_mn(2,1);

c_b = sqrt(1.8*f*h*c_ql);

% % figure()
% % plot(f,c_b,'LineWidth',1.2);
% % grid minor
% % xlabel('f [Hz]','Fontsize',15); 
% % ylabel('c_b [m/s]','Fontsize',15);

%% e) first five modal frequencies

f_b_mn = 0.4694*c_ql*h/r_m^2*[1 2.08 3.41 3.89 5];

%% DATA ex3

rho_s = 5000; %[kg/m^3]
r_s = 0.001; %[m]
L = 0.4; %[m]

V_s = L*pi*r_s^2;
mass_s = rho_s*V_s;
S = pi*r_s^2;

Q_p = 50;

E_s = 200e9; % [Pa]
%% g) tension to tune string and soundboard

f_s = f_b_mn(1);

%f_s = 1/(2*L)*sqrt(T/(rho_s*pi*r_s^2));
T = rho_s*pi*r_s^2*4*L^2*f_s^2;

%% h) freqs of first five modes considering stiffness

gyration_radius = r_s/2;
%inertia_moment = m*R^2/2;

B = pi * E_s * S * gyration_radius^2/(T*L^2);

mu = mass_s/L;
c_s = sqrt(T/mu);

n = [1 2 3 4 5];
f_n = n.*c_s/(2*L);
f_n_stiffness = n .* f_n(1) .* sqrt(1 + (B.*n.^2) )*(1 + 2/pi * B^(0.5) + (4/pi^2)*B);

%% i) freqs of modes considering coupling string-plate (neglecting stiffness)

m = linspace(1,5,5);
coupling_coeff = mass_s./(m.^2*mass_p)*10^3;

w_n = 2*pi*f_n;
w_b_mn = 2*pi*f_b_mn;
diff = (w_n - w_b_mn)./(w_b_mn);

omega_plus = zeros(1,5);
omega_minus = zeros(1,5);
vect = [0.077 -0.045 0.012 -0.13 0.001 -0.015 0.04 -0.02 0.02];

omega_plus(1) = vect(1)*w_n(1)/2;
omega_minus(1) = omega_plus(1);

for i=1:4
    omega_minus(i+1) = vect(2*i)*w_b_mn(i+1) + 1;
    omega_plus(i+1) = vect(2*i+1)*w_b_mn(i+1) + 1;
end

