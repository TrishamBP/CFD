% PSOR and LSOR for Laplace's PDE in 2d  with fixed BC

clear all
%1 is for PSOR (omega=1 becomes PGS)
%2 is for  LSOR 

method =2;

maxiter=1000;
convergence=1.e-5;

Length=1;
Height=2;
dx=0.05;
dy= 0.05;
beta=dx/dy;
nx=Length/dx +1;
ny=Height/dy + 1;

for i=1:nx
for j=1:ny
    x(i,j)=dx*(i-1);
    y(i,j)=dy*(j-1);
end
end


if method==1
for kom=1:6
omega(kom)=kom*0.3 

%Implement Bounday Conditions

u_old(1:nx,1)=100;
u_old(1:nx,ny)=0;
u_old(1,1:ny)=0;
u_old(nx,1:ny)=0;
u_new(1:nx,1)=100;
u_new(1:nx,ny)=0;
u_new(1,1:ny)=0;
u_new(nx,1:ny)=0;

%Initial guess
for i=2:nx-1
    for j=2:ny-1
        u_old(i,j)=0;
        u_new(i,j)=0;
        
    end
end

n=1;
sum_error=1.e6;
iter(kom,1)=1;
  while sum_error>convergence
    sum_error=0;
    for j = 2:(ny-1)
        for i=2:(nx-1)
            xx=u_old(i+1,j)+u_new(i-1,j)+beta^2*u_old(i,j+1)+beta^2*u_new(i,j-1)-2*(1+beta^2)*u_old(i,j);     
            u_new(i,j)=u_old(i,j) +(omega(kom)/(2*(1+beta^2))) *xx;
            sum_error = sum_error + abs(u_new(i,j)-u_old(i,j));
            u_old(i,j)=u_new(i,j);
        end
    end
    error(kom,n)=sum_error;
    iter(kom,n)=n;
    n=n+1;
  end
end
save('runhmw5','omega')
appendedstring1=['\omega= ',num2str(omega(1)),''];
appendedstring2=['\omega= ',num2str(omega(2)),''];
appendedstring3=['\omega= ',num2str(omega(3)),''];
appendedstring4=['\omega= ',num2str(omega(4)),''];
appendedstring5=['\omega= ',num2str(omega(5)),''];
appendedstring6=['\omega= ',num2str(omega(6)),''];

subplot(2,2,1),semilogy(iter(1,:),error(1,:),iter(2,:),error(2,:),iter(3,:),error(3,:),iter(4,:),error(4,:),iter(5,:),error(5,:),iter(6,:),error(6,:))
legend(appendedstring1,appendedstring2,appendedstring3,appendedstring4,appendedstring5,appendedstring6)
title('Convergence History')
xlabel('Number of Iterations')
ylabel('Error')

%Plot the solution with the last Omega

subplot(2,2,2),[cs,h]=contour(x,y,u_new);
clabel(cs,h);
title('Numerical Solution')
xlabel('X(m)')
ylabel('Y(m)')

%calculates the analytical solution and the solution error using the last numerical solution
for i=1:nx
for j=1:ny
    u_ana(i,j)=0;
    for n=1:20
        u_ana(i,j)= u_ana(i,j)+(1-(-1)^n)*(sinh(n*pi*(Height-y(i,j))/Length)/sinh(n*pi*Height/Length))*(sin(n*pi*x(i,j)/Length))/(n*pi);
    end
    u_ana(i,j)=100*2*u_ana(i,j);
 sol_error(i,j)=abs(u_ana(i,j)-u_new(i,j));
end
end

subplot(2,2,3),[cs,h]=contour(x,y,u_ana);
clabel(cs,h);
title('Analytical Solution')
xlabel('X(m)')
ylabel('Y(m)')

z(1)=1.;
z(2)=0.1;
z(3)=0.01;
subplot(2,2,4),[cs,h]=contour(x,y,sol_error,z);
clabel(cs,h);
title('Error')

end

if method==2
    
  for kom=1:7
    omega(kom)=1.18+kom*0.02

%Implement Boundary Conditions

u_old(1:nx,1)=100;
u_old(1:nx,ny)=0;
u_old(1,1:ny)=0;
u_old(nx,1:ny)=0;
u_new(1:nx,1)=100;
u_new(1:nx,ny)=0;
u_new(1,1:ny)=0;
u_new(nx,1:ny)=0;

%Initial guess
for i=2:nx-1
    for j=2:ny-1
        u_old(i,j)=0;
        u_new(i,j)=0;
        
    end
end

n=1;
sum_error=1.e6;
iter(kom,1)=1;
  while sum_error>convergence  
      sum_error=0;

   for j=2:ny-1
       %form the Tridiagonal matrix for the nx-2 inner grid points
       at=omega(kom);
      bt=-2*(1+beta^2);
      ct=omega(kom);
      e=ones(nx-2,1);
      T=spdiags([at*e,bt*e,ct*e],-1:1,nx-2,nx-2);
      %fill the rhs of the equation
      rhs(1,1)=-omega(kom)*u_new(1,j)-(1-omega(kom))*(2+2*beta^2)*u_old(2,j)-omega(kom)*beta^2*(u_old(2,j+1)+u_new(2,j-1));
      rhs(nx-2,1)=-omega(kom)*u_new(nx,j)-(1-omega(kom))*(2+2*beta^2)*u_old(nx-1,j)-omega(kom)*beta^2*(u_old(nx-1,j+1)+u_new(nx-1,j-1));
      
      for ii=2:nx-3
          rhs(ii,1)=-(1-omega(kom))*2*(1+beta^2)*u_old(ii+1,j)-omega(kom)*beta^2*(u_old(ii+1,j+1)+u_new(ii+1,j-1));
      end
      [L,U,P]=lu(T);
      sol=U\(L\(P'*rhs));
      for ii=1:nx-2
        u_new(ii+1,j)=sol(ii);
        sum_error = sum_error + abs(u_new(ii+1,j)-u_old(ii+1,j));
        u_old(ii+1,j)=u_new(ii+1,j);
      end
   end
   
   error(kom,n)=sum_error;
    iter(kom,n)=n;
    n=n+1;
  end
  end


appendedstring1=['\omega= ',num2str(omega(1)),''];
appendedstring2=['\omega= ',num2str(omega(2)),''];
appendedstring3=['\omega= ',num2str(omega(3)),''];
appendedstring4=['\omega= ',num2str(omega(4)),''];
appendedstring5=['\omega= ',num2str(omega(5)),''];
appendedstring6=['\omega= ',num2str(omega(6)),''];
appendedstring7=['\omega= ',num2str(omega(7)),''];

figure(1),semilogy(iter(1,:),error(1,:),iter(2,:),error(2,:),iter(3,:),error(3,:),iter(4,:),error(4,:),iter(5,:),error(5,:),iter(6,:),error(6,:),iter(7,:),error(7,:));
legend(appendedstring1,appendedstring2,appendedstring3,appendedstring4,appendedstring5,appendedstring6,appendedstring7);

%plot the solution. Note that u_new is obtained with Omega(7)

figure(2),[cs,h]=contour(x,y,u_new);
clabel(cs,h);

%calculates the analytical solution and the solution error using the last
%numerical solution as u_new

for i=1:nx
for j=1:ny
    u_ana(i,j)=0;
    for n=1:20
       u_ana(i,j)= u_ana(i,j)+(1-(-1)^n)*(sinh(n*pi*(Height-y(i,j)/Length))/sinh(n*pi*Height/Length))*(sin(n*pi*x(i,j)/Length))/(n*pi);
    end
    u_ana(i,j)=100*2*u_ana(i,j);
 sol_error(i,j)=abs(u_ana(i,j)-u_new(i,j));
end
end

figure(3),[cs,h]=contour(x,y,u_ana);
clabel(cs,h);

z(1)=1.;
z(2)=0.1;
z(3)=0.01;
figure(4),[cs,h]=contour(x,y,sol_error,z);
clabel(cs,h);

end


