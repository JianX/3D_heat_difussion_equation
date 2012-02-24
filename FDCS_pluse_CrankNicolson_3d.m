function FDCS_pluse_CrankNicolson_3d(nx)
%%
% I made a modification for the first time.
% 1. this program solve the heat diffusion equation in 3 dimensions. First, it is set up on an
% arbitrary Cartesian domain. And dx != dy != dz. 
% 2. It implements the following two methods of solving the system:
% - Forward-Time Central-Space 
% - Crank-Nicolson Method (20 pts) 
% 3. My code must have the following features built in:
% - A function that sets the initial condition (Gaussian+noise) 
% - Ability to specify the following boundary conditions. 
% - Dirichlet (constant) conditions
% - An arbitrary time-independent "source term" added to the right hand side of the equation
% (ie specified S(i,j,k), where S is given)

%%
%set up time step and grid

nsteps=100000;
alpha=0.0005;
dt=0.02;
dx=0.04;
dy=0.05;
dz=0.04;
c1=alpha*dt/dx^2;
c2=alpha*dt/dy^2;
c3=alpha*dt/dz^2;

xx=linspace(0,nx*dx,nx);
yy=linspace(0,nx*dy,nx);
zz=linspace(0,nx*dz,nx);
[x y z]=meshgrid(xx,yy,zz); %create grid
%%
%initial contidition with Guassian and noise
T=exp(-(x-nx*dx/2).^2).*exp(-(y-nx*dy/2).^2).*exp(-(z-nx*dz/2).^2); %create gaussian profile
v=ones(nx,nx,nx);
R = (2*normrnd(0,v)-1)/10;% Gaussian noise for 3D
T=T+R;
%add sourse term
S = sin(x).*cos(z).*tan(y);
Sreshape=reshape(S,1,nx^3)';
%initial plot
% subplot(1,2,1)
% imagesc(T(:,:,round(nx/2)))
% set(gca,'ZLim',[0,1]);
% subplot(1,2,2)
% imagesc(T(:,:,round(nx/3)))
% set(gca,'ZLim',[0,1]);

%%
%set up Dirichlet boundary condition
X1=0;XL=0; %this boundary condition could be any constant
Y1=0;YL=0;
Z1=0;ZL=0;
T(1,:,:)=X1;T(nx,:,:)=XL;
T(:,1,:)=Y1;T(:,nx,:)=YL;
T(:,:,1)=Z1;T(:,:,nx)=ZL;
%%
%find the matrix for FDCS
A=zeros(nx,1);
B=zeros(nx);
C=zeros(nx,nx,nx);
for k=2:nx-1
    A(k)=-2*c1-2*c2-2*c3;
end
for j=2:nx-1
    B(j,:)=A;
end
for i=2:nx-1
    C(:,:,i)=B;
end
media=C(:);
leng=size(media,1);
C1=zeros(leng,1);
C2=C1;
C3=C1;
C11=C1;
C22=C1;
C33=C1;
%for FDCS, produce vectors for spdiags
C1(nx^2+nx+1:nx^3-nx^2-nx-2)=c1;
C2(nx^2+2:nx^3-nx^2-2*nx-1)=c2;
C3(nx+2:nx^3-2*nx^2-nx-1)=c3;

C11(nx^2+nx+1:nx^3-nx^2-nx-2)=c1;
C22(nx^2+2:nx^3-nx^2-2*nx-1)=c2;
C33(nx+2:nx^3-2*nx^2-nx-1)=c3;
C11=flipud(C11);
C22=flipud(C22);
C33=flipud(C33);

%put all the vectors in a matrix B for spdiags
B=zeros(nx^3);
B(:,1)=C3;
B(:,2)=C2;
B(:,3)=C1;
B(:,4)=C(:);
B(:,5)=C11;
B(:,6)=C22;
B(:,7)=C33;
Dvector=[-(nx^2) -nx -1 0 1 nx nx^2];
Eureka = spdiags(B,Dvector,nx^3,nx^3);
%Eureka is the kernel for three matrix
Eureka_CN1=Eureka/2+diag(ones(nx^3,1),0);
Eureka_CN2=-Eureka/2+diag(ones(nx^3,1),0);
Eureka=Eureka+diag(ones(nx^3,1),0);

save Eureka;
save Eureka_CN1;
save Eureka_CN2;

%%
%solve it
for i=1:nsteps
    b=reshape(T,1,nx^3)';
    %FDCS method
    b=Eureka*b+dt*Sreshape;
    b_CN=Eureka_CN2\(Eureka_CN1*b+dt*Sreshape);
    T=reshape(b,nx,nx,nx);
    T_CN=reshape(b_CN,nx,nx,nx);
    T(1,:,:)=X1;T(nx,:,:)=XL;
    T(:,1,:)=Y1;T(:,nx,:)=YL;
    T(:,:,1)=Z1;T(:,:,nx)=ZL;
    T_CN(1,:,:)=X1;T_CN(nx,:,:)=XL;
    T_CN(:,1,:)=Y1;T_CN(:,nx,:)=YL;
    T_CN(:,:,1)=Z1;T_CN(:,:,nx)=ZL;
    if (~mod(i,1))
        subplot(1,2,1)
        imagesc(T(:,:,round(nx/2)))
        set(gca,'ZLim',[0,1]);
        title('FDCS')
        subplot(1,2,2)
        imagesc(T_CN(:,:,round(nx/2)))
        set(gca,'ZLim',[0,1]);
        title('CrankNicolson')
        pause(.001);
    end
end

















