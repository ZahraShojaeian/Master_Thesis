n=3/2;
phi_prime= pi/2;
beta_prime= 2*pi/9;
 e_diff= zeros(1,m);
 E_diff_xx2= zeros(1,m);
for s=1:1*m
    for r=1:1*m   
    for j=1:m
        if x(s)>x_prime(j)
            beta(j)=atan(h/(x(s)-x_prime(j))) + pi/2;
        elseif x_prime(j)>=x(s)
            beta(j)=atan((x_prime(j)-x(s))/h);
        end
        if y(r)>y_prime(1)
            phi(j)=atan(h/(-y_prime(1)+y(r)));
            
        elseif y_prime(1)>=y(r)
             phi(j)=atan((-y(r)+y_prime(1))/h) + pi/2;
        end
         alpha1(j)= asin((sqrt(sin(beta_prime)^2-((sin(beta(j))^2)*(cos(phi(j))^2))))/sin(beta_prime));
         alpha2(j)= asin((sqrt(sin(beta_prime)^2-((sin(beta(j))^2)*(cos((n*pi)-phi(j))^2))))/sin(beta_prime)); 
         Q= 2*sin(0.5*(beta(j)-beta_prime))*sin(0.5*(beta(j)+beta_prime))*(1+(cos(beta(j))*cos(beta_prime)))/(sin(beta(j))*sin(beta_prime));
         De(j)= ((1/n)*(sin(phi_prime/n)))/(cos((pi-alpha1(j))/n)-cos(phi_prime/n)) + ((1/n)*(sin(phi_prime/n)))/(cos((pi-alpha2(j))/n)+cos(phi_prime/n));
         De_prime= -(sin(phi_prime))/(cos(alpha1(j))+cos(phi_prime)); 
         Ie(j) = Ex(j,m)*2*(De(j)-De_prime)/(1i*k*etha*(sin(beta_prime)^2));

         e_diff(j)= etha*(1)*Ie(j)*((exp(-1*i*k*(sqrt((((x(s))-x_prime(j))^2)+((y(r)-y_prime(1))^2)+((h)^2)))))/(4*pi*sqrt((((x(s))-x_prime(j))^2)+((y(r)-y_prime(1))^2)+((h)^2))));
         if  beta(j)> (pi-beta_prime)
    e_diff(j)=0;
end
    end
    E_diff_xx2(r,s)= 1i*k*sum(e_diff);
    end
end
figure
surf(xx,yy,zz,10*log10(abs(Ex)), 'EdgeColor', 'none' )
 view(2)
figure
surf(x,y,z,10*log10(0.05*abs(E_diff_xx1+E_diff_xx2)), 'EdgeColor', 'none' )
 view(2)
 figure
surf(x,y,z,10*log10(0.05*abs(E_diff_xx1-E_diff_xx2)), 'EdgeColor', 'none' )
 view(2)
 figure
surf(x,y,z,10*log10(0.05*abs(E_diff_xx1)), 'EdgeColor', 'none' )
 view(2)
 figure
surf(x,y,z,10*log10(0.05*abs(E_diff_xx2)), 'EdgeColor', 'none' )
 view(2)

%%%%%% Image reconstruction using Quadratic phase approximation and Fresnel Diffraction Integral %%%%%%

EX= (E_diff_xx1-E_diff_xx2);
X=linspace(-0.15,0.15,1*m);
Y=linspace(-0.15,0.15,1*m);
XX=linspace(-0.15,0.15,1*m);
YY=linspace(-0.3,0,1*m);
Image_2= zeros(1,m);
Z=0.4;
p=(1/Z)-(1/0.3);
for u=1:1*m
    for v=1:1*m
        for s=1:1*m
            for r=1:1*m
        if ((X(s)^2)+(Y(r)^2))> 0.0225
            EX(s,r)=0;
        end
            Image_1(s,r)= EX(s,r)*exp(1i*(X(s)^2 + Y(r)^2)*p/2)*exp(-1i*2*pi*((X(s)*XX(u))+(Y(r)*YY(v)))/(lambda*Z));
        
            end
        end
        Image(u,v)=( exp(((XX(u)^2) + (YY(v)^2))*1i*k/(2*Z))/(1i*lambda*Z)).*(sum(Image_1(:)));
    end
end
figure
surf(XX,YY,z,10*log10((abs(Image)/abs(max(max(Image))))), 'EdgeColor', 'none' )
 view(2)
  figure
surf(XX,YY,z,10*log10((abs(Image))), 'EdgeColor', 'none' )
 view(2)