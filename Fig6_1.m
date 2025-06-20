clear 
clf
Smax=200;Smin=0;Vmax=1;Vmin=0;I=100;J=100;dt=1/12000;  
r=0.05;kappa=2.5;theta=0.16;rho=0.1;sigma=0.45;T=1;K=100;    
N=1+ceil(T/dt);% dt = T/nt;  B = ceil(A) rounds the elements of A to the nearest integers greater than or equal to A. 
ds = (Smax-Smin)/(I); % step length of s 
dv = (Vmax-Vmin)/(J); % step length of v  
U=zeros(I+1,J+1,N); 
for i=1:I+1   
for j=1:J+1
U(i,j,1)=max(0,(i-1)*ds-K);  
end
end

k=0; deltat=1/6;
for n=1:N-1 % the interior elements. Cross term part.
for j=2:J   
for i = 2:I

U(i,j,n+1)=U(i,j,n)+((j-1)*dv)*((i-1)*ds)^2*dt/(2*ds^2)*(U(i+1,j,n)-2*U(i,j,n)+U(i-1,j,n))...
+sigma^2*((j-1)*dv)^3*dt/(2*dv^2)*(U(i,j+1,n)-2*U(i,j,n)+U(i,j-1,n))...
+r*(i-1)*ds*dt/(2*ds)*(U(i+1,j,n)-U(i-1,j,n))...
+kappa*(theta-(j-1)*dv)*(j-1)*dv*dt/(2*dv)*(U(i,j+1,n)-U(i,j-1,n))...
+rho*sigma*((j-1)*dv)^2*((i-1)*ds)*dt/(4*ds*dv)*(U(i+1,j+1,n)+U(i-1,j-1,n)-U(i-1,j+1,n)-U(i+1,j-1,n))...
-r*dt*U(i,j,n);
end
end
for j=1:J+1
U(1,j,n+1)=0;  %%  %%S=0 U(0,sigma,t)=0  
end
for i=2:I
U(i,1,n+1)= r*(i-1)*dt*(U(i+1,1,n)-U(i,1,n))+U(i,1,n)*(1-r*dt); %%V=0
end
for i=2:I
U(i,J+1,n+1)=(i-1)*ds; %%V=Vmax
end
for j=1:J+1
U(I+1,j,n+1)=ds+U(I,j,n+1);  %% %%\partial U(S_max,sigma,t)/(\partial S)  =1
end
end
U0=U;
WTS15=U0(41,:,n);%%% S=80 V=0:1
WTV25=U0(:,26,n);%%% 
k=0.005; deltat=1/6;
for n=1:N-1 % the interior elements. Cross term part.
for j=2:J   
for i = 2:I
f1=(U(i+1,j,n)-2*U(i,j,n)+U(i-1,j,n))/(ds)^2; %% S^2
f2=(U(i+1,j+1,n)+U(i-1,j-1,n)-U(i-1,j+1,n)-U(i+1,j-1,n))/(4*ds*dv);  %% SV
f3=(U(i,j+1,n)-2*U(i,j,n)+U(i,j-1,n))/(dv)^2; %% V^2
F1=k*((i-1)*ds)*dt*sqrt(2/(pi*deltat))*sqrt(((j-1)*dv)*((i-1)*ds)^2*(f1)^2 +...
sigma^2*((j-1)*dv)^3*(f2)^2+2*rho*sigma*((j-1)*dv)^2*((i-1)*ds)*(f2*f1));%%cost1

U(i,j,n+1)=U(i,j,n)+((j-1)*dv)*((i-1)*ds)^2*dt/(2*ds^2)*(U(i+1,j,n)-2*U(i,j,n)+U(i-1,j,n))...
+sigma^2*((j-1)*dv)^3*dt/(2*dv^2)*(U(i,j+1,n)-2*U(i,j,n)+U(i,j-1,n))...
+r*(i-1)*ds*dt/(2*ds)*(U(i+1,j,n)-U(i-1,j,n))...
+kappa*(theta-(j-1)*dv)*(j-1)*dv*dt/(2*dv)*(U(i,j+1,n)-U(i,j-1,n))...
+rho*sigma*((j-1)*dv)^2*((i-1)*ds)*dt/(4*ds*dv)*(U(i+1,j+1,n)+U(i-1,j-1,n)-U(i-1,j+1,n)-U(i+1,j-1,n))...
-r*dt*U(i,j,n)-F1;
end
end
for j=1:J+1
U(1,j,n+1)=0;  %%  %%S=0 U(0,sigma,t)=0   i=1=S0
end
for i=2:I
U(i,1,n+1)= r*(i-1)*dt*(U(i+1,1,n)-U(i,1,n))+U(i,1,n)*(1-r*dt); %%V=0
end
for i=2:I
U(i,J+1,n+1)=(i-1)*ds; %%V=Vmax U(S,Vmax,t)=S
end
for j=1:J+1
U(I+1,j,n+1)=ds+U(I,j,n+1);  %% %%\partial U(S_max,sigma,t)/(\partial S)  =1 
end
end
U1=U;
WTS15K001=U1(41,:,n);%%% S=80 V=0:1
WTV25K001=U1(:,26,n);%%%   V=0.1,S=0:50,k=0.1

k=0.01; deltat=1/6;
for n=1:N-1 % the interior elements. Cross term part.
for j=2:J   %%
for i = 2:I
f1=(U(i+1,j,n)-2*U(i,j,n)+U(i-1,j,n))/(ds)^2; %% S^2
f2=(U(i+1,j+1,n)+U(i-1,j-1,n)-U(i-1,j+1,n)-U(i+1,j-1,n))/(4*ds*dv);  %% SV
f3=(U(i,j+1,n)-2*U(i,j,n)+U(i,j-1,n))/(dv)^2; %% V^2
F1=k*((i-1)*ds)*dt*sqrt(2/(pi*deltat))*sqrt(((j-1)*dv)*((i-1)*ds)^2*(f1)^2 +...
sigma^2*((j-1)*dv)^3*(f2)^2+2*rho*sigma*((j-1)*dv)^2*((i-1)*ds)*(f2*f1));%%cost1

U(i,j,n+1)=U(i,j,n)+((j-1)*dv)*((i-1)*ds)^2*dt/(2*ds^2)*(U(i+1,j,n)-2*U(i,j,n)+U(i-1,j,n))...
+sigma^2*((j-1)*dv)^3*dt/(2*dv^2)*(U(i,j+1,n)-2*U(i,j,n)+U(i,j-1,n))...
+r*(i-1)*ds*dt/(2*ds)*(U(i+1,j,n)-U(i-1,j,n))...
+kappa*(theta-(j-1)*dv)*(j-1)*dv*dt/(2*dv)*(U(i,j+1,n)-U(i,j-1,n))...
+rho*sigma*((j-1)*dv)^2*((i-1)*ds)*dt/(4*ds*dv)*(U(i+1,j+1,n)+U(i-1,j-1,n)-U(i-1,j+1,n)-U(i+1,j-1,n))...
-r*dt*U(i,j,n)-F1;
end
end
for j=1:J+1
U(1,j,n+1)=0;  %%  %%S=0 U(0,sigma,t)=0   
end
for i=2:I
U(i,1,n+1)= r*(i-1)*dt*(U(i+1,1,n)-U(i,1,n))+U(i,1,n)*(1-r*dt); %%V=0
end
for i=2:I
U(i,J+1,n+1)=(i-1)*ds; %%V=Vmax
end
for j=1:J+1
U(I+1,j,n+1)=ds+U(I,j,n+1); 
end
end
U2=U;
WTS15K002=U2(41,:,n);%%% S=80 V=0:1
WTV25K002=U2(:,26,n);

k=0.015; deltat=1/6;
for n=1:N-1 % the interior elements. Cross term part.
for j=2:J   
for i = 2:I
f1=(U(i+1,j,n)-2*U(i,j,n)+U(i-1,j,n))/(ds)^2; %% S^2
f2=(U(i+1,j+1,n)+U(i-1,j-1,n)-U(i-1,j+1,n)-U(i+1,j-1,n))/(4*ds*dv);  %% SV
f3=(U(i,j+1,n)-2*U(i,j,n)+U(i,j-1,n))/(dv)^2; %% V^2
F1=k*((i-1)*ds)*dt*sqrt(2/(pi*deltat))*sqrt(((j-1)*dv)*((i-1)*ds)^2*(f1)^2 +...
sigma^2*((j-1)*dv)^3*(f2)^2+2*rho*sigma*((j-1)*dv)^2*((i-1)*ds)*(f2*f1));%%cost1

U(i,j,n+1)=U(i,j,n)+((j-1)*dv)*((i-1)*ds)^2*dt/(2*ds^2)*(U(i+1,j,n)-2*U(i,j,n)+U(i-1,j,n))...
+sigma^2*((j-1)*dv)^3*dt/(2*dv^2)*(U(i,j+1,n)-2*U(i,j,n)+U(i,j-1,n))...
+r*(i-1)*ds*dt/(2*ds)*(U(i+1,j,n)-U(i-1,j,n))...
+kappa*(theta-(j-1)*dv)*(j-1)*dv*dt/(2*dv)*(U(i,j+1,n)-U(i,j-1,n))...
+rho*sigma*((j-1)*dv)^2*((i-1)*ds)*dt/(4*ds*dv)*(U(i+1,j+1,n)+U(i-1,j-1,n)-U(i-1,j+1,n)-U(i+1,j-1,n))...
-r*dt*U(i,j,n)-F1;
end
end
for j=1:J+1
U(1,j,n+1)=0;  %%  %%S=0 U(0,sigma,t)=0   
end
for i=2:I
U(i,1,n+1)= r*(i-1)*dt*(U(i+1,1,n)-U(i,1,n))+U(i,1,n)*(1-r*dt); %%V=0
end
for i=2:I
U(i,J+1,n+1)=(i-1)*ds;
end
for j=1:J+1
U(I+1,j,n+1)=ds+U(I,j,n+1);  %% %%\partial U(S_max,sigma,t)/(\partial S)  =1 
end
end
U3=U;

k=0.02; deltat=1/6;
for n=1:N-1 % the interior elements. Cross term part.
for j=2:J   
for i = 2:I
f1=(U(i+1,j,n)-2*U(i,j,n)+U(i-1,j,n))/(ds)^2; %% S^2
f2=(U(i+1,j+1,n)+U(i-1,j-1,n)-U(i-1,j+1,n)-U(i+1,j-1,n))/(4*ds*dv);  %% SV
f3=(U(i,j+1,n)-2*U(i,j,n)+U(i,j-1,n))/(dv)^2; %% V^2
F1=k*((i-1)*ds)*dt*sqrt(2/(pi*deltat))*sqrt(((j-1)*dv)*((i-1)*ds)^2*(f1)^2 +...
sigma^2*((j-1)*dv)^3*(f2)^2+2*rho*sigma*((j-1)*dv)^2*((i-1)*ds)*(f2*f1));%%cost1

U(i,j,n+1)=U(i,j,n)+((j-1)*dv)*((i-1)*ds)^2*dt/(2*ds^2)*(U(i+1,j,n)-2*U(i,j,n)+U(i-1,j,n))...
+sigma^2*((j-1)*dv)^3*dt/(2*dv^2)*(U(i,j+1,n)-2*U(i,j,n)+U(i,j-1,n))...
+r*(i-1)*ds*dt/(2*ds)*(U(i+1,j,n)-U(i-1,j,n))...
+kappa*(theta-(j-1)*dv)*(j-1)*dv*dt/(2*dv)*(U(i,j+1,n)-U(i,j-1,n))...
+rho*sigma*((j-1)*dv)^2*((i-1)*ds)*dt/(4*ds*dv)*(U(i+1,j+1,n)+U(i-1,j-1,n)-U(i-1,j+1,n)-U(i+1,j-1,n))...
-r*dt*U(i,j,n)-F1;
end
end
for j=1:J+1
U(1,j,n+1)=0;  %%  %%S=0 U(0,sigma,t)=0   
end
for i=2:I
U(i,1,n+1)= r*(i-1)*dt*(U(i+1,1,n)-U(i,1,n))+U(i,1,n)*(1-r*dt); 
end
for i=2:I
U(i,J+1,n+1)=(i-1)*ds; %%
end
for j=1:J+1
U(I+1,j,n+1)=ds+U(I,j,n+1);
end
end
U4=U;

S=Smin:ds:Smax;
V=Vmin:dv:Vmax;
x=S;
y=V;
[xi,yj]=meshgrid(x,y);

subplot(2,2,1)
mesh(S,V,U0(:,:,N)-U1(:,:,N));  % n=N-1  U(i,j,n+1)
pcolor(xi,yj,U0(:,:,N)-U1(:,:,N));
shading interp;    %%
colorbar; colormap(jet);
xlabel('S','fontsize',16,'FontWeight','bold')
ylabel('V','fontsize',16,'FontWeight','bold')
title('Difference in option prices, when k_{TC}=0.005','fontsize',16,'FontWeight','bold')

subplot(2,2,2)
mesh(S,V,U0(:,:,N)-U2(:,:,N));  % n=N-1  U(i,j,n+1)
pcolor(xi,yj,U0(:,:,N)-U2(:,:,N));
shading interp;    %%
colorbar; colormap(jet);
xlabel('S','fontsize',16,'FontWeight','bold')
ylabel('V','fontsize',16,'FontWeight','bold')
title('Difference in option prices, when k_{TC}=0.01','fontsize',16,'FontWeight','bold')

subplot(2,2,3)
mesh(S,V,U0(:,:,N)-U3(:,:,N));  % n=N-1  U(i,j,n+1)  
pcolor(xi,yj,U0(:,:,N)-U3(:,:,N));
shading interp;    
colorbar; colormap(jet);
xlabel('S','fontsize',16,'FontWeight','bold')
ylabel('V','fontsize',16,'FontWeight','bold')
title('Difference in option prices, when k_{TC}=0.015','fontsize',16,'FontWeight','bold')
subplot(2,2,4)
mesh(S,V,U0(:,:,N)-U4(:,:,N));  % n=N-1  U(i,j,n+1)
pcolor(xi,yj,U0(:,:,N)-U4(:,:,N));
shading interp;    
colorbar; colormap(jet);
xlabel('S','fontsize',16,'FontWeight','bold')
ylabel('V','fontsize',16,'FontWeight','bold')
title('Difference in option prices, when k_{TC}=0.02','fontsize',16,'FontWeight','bold')
