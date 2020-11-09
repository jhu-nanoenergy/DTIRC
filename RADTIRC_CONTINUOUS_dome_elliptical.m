n=1.41; %refractive index
K=0.3; % a^2/b^2 of the ellipse.
phi=-pi/39.2; %predefined arc angle at x-axis
thetaa1=pi/7.2; % acceptance angle at x-axis
thetaa2=pi/8; % acceptance angle at y-axis
d00=5; %size of the bottom square aperture
ksi=atan2(K*sin(pi/2+phi),cos(pi/2+phi))-pi/2; % vertical incidental angle
Theta=asin(sin(thetaa1-ksi)/n)+ksi; % maximum vertical ray angle
thetac=asin(1/n); %critical angle 
theta0=pi-Theta-2*thetac; %exiting angle of parallel rays
d1=n*sin(theta0)*d00/sin(thetaa1); %entrance aperture size by etuende conservation
b=d1/2*sqrt((tan(phi)^2+K)/(K*tan(phi)^2)); % semi-major axis
a=b*sqrt(K);
h=sign(phi)*b*(1-sqrt(K/(tan(phi)^2+K))); %height of the dome
H0=0.5*(d1+d00)*cot(Theta)+h; %height of the collector.
PC=[];
S=[];
dome_area=0;
thetaa=thetaa1;


x_3D=[];
y_3D=[];
z_3D=[];
for range=0:7
    for alpha=range*pi/4:pi/4/8:(range+1)*pi/4 %rotation angle
        %     thetaa1_temp=thetaa1;
        %
        %
        %     thetaa2_temp=asin(tan(alpha)*sin(thetaa1));
        %     if(thetaa2_temp>thetaa2)
        %         thetaa2_temp=thetaa2;
        %         thetaa1_temp=asin(tan(pi/2-alpha)*sin(thetaa2));
        %     end
        %     thetaa=2*asin(sqrt((2*sin(thetaa1_temp/2))^2+(2*sin(thetaa2_temp/2))^2)/2);
        
%         
%         thetaa=thetaa1;
        id=(mod(ceil(range/2),2)==0);
        d0=abs(d00/(cos(alpha)*id+sin(alpha)*(~id))); % exit aperture
        
%             thetaa=thetaa2;
     
        
        fun=@(x) funfun(x,n,K,b,d0,H0); %solve for the arc angle
        sols=fsolve(fun,[phi,thetaa]);
        phi=sols(1);
        thetaa=sols(2);
        
        ksiM=atan2(K*sin(pi/2+phi),cos(pi/2+phi))-pi/2;
        theta=linspace(-phi,phi,160);
        ksi=atan2(K*sin(pi/2+theta),cos(pi/2+theta))-pi/2;
        thetap=asin(sin(thetaa-ksi)/n)+ksi;
        Theta=asin(sin(thetaa-ksiM)/n)+ksiM;  % maximum vertical ray angle
        thetac=asin(1/n);
        theta0=pi-Theta-2*thetac;
        
        d1=n*sin(theta0)*d0/sin(thetaa);
        
        x_phi=-b*sqrt(K*tan(phi)^2/(K+tan(phi)^2));
        z_phi=sign(phi)*b*sqrt(K/(K+tan(phi)^2));
        x_theta=-b*sign(phi)*tan(theta).*sqrt(K./(K+tan(theta).^2));
        z_theta=sign(phi)*b*sqrt(K./(K+tan(theta).^2));
        dist=sqrt((x_theta-x_phi).^2+(z_theta-z_phi).^2); % distance to the incident points of max ray
        l1=dist.*sin(thetaa-asin((z_theta-z_phi)./dist));   
        l1(end)=0;
        %segment 1 of the optical path
        h=sign(phi)*b*(1-sqrt(K/(tan(phi)^2+K)));
        H=0.5*(d1+d0)*cot(Theta)+h;
        
        C=2*abs(x_phi)*sin(thetaa)+n*(d1+d0)/(2*sin(Theta)); % constant optical path
   
        
        
        y0=sign(phi)*b-z_theta;
        U=(C-l1)/n;
        V=x_theta+d0/2;
        W=H-y0;
        l2=(U.^2-V.^2-W.^2)./(2*(U+V.*sin(thetap)-W.*cos(thetap))); %segment 2 for rays hitting the same point
        x=x_theta+l2.*sin(thetap);
        y=H-l2.*cos(thetap)-y0;
              
        
        l2=(U-V*sin(theta0)-W*cos(theta0))./(2*sin((theta0+thetap)/2).^2); %segment 2 for rays exiting parallelly 
        xp=x_theta+l2.*sin(thetap);
        yp=H-l2.*cos(thetap)-y0;
        
                      
        [X,XP]=meshgrid(x,xp); %findin the intersection of 2 choices of segment 2
        [Y,YP]=meshgrid(y,yp);
        D=(X-XP).^2+(Y-YP).^2;
        Dmin=min(min(D));
        [~,t]=find(Dmin==D); %find the intersection of two side segments
        y0=y(t); % the intersection 
        x_dome=linspace(0,d1/2,20);
        y_dome=b*(sqrt(1-x_dome.^2/a^2)-1)+H;
        xb=linspace(0.02,d0/2,20); %coordinates for exit aperture
        yb=zeros(1,20);
        xf=[(y>y0).*x+(y<=y0).*xp,x_dome,xb]; %The side profile at a specific angle, including the dome
        yf=[(y>y0).*y+(y<=y0).*yp,y_dome,yb];
       
%         if alpha<eps||(abs(alpha-pi/4)<eps)||(abs(alpha-pi/2)<eps)
%             %generate the bracket plane at 0, pi/4
%             for k=1:length(yf)
%                 x_temp=linspace(0,d1,20);
%                 x_side=x_temp(x_temp<=xf(k));
%                 xf=[xf,x_side];
%                 yf=[yf,repmat(yf(k),1,length(x_side))];
%             end
%         end
        
        
        x_new=xf*cos(alpha); %Coordinates Transformation
        z_new=xf*sin(alpha);
        
       
     
        
        x_3D=[x_3D,x_new]; %Append the entire point sets with the side profile
        y_3D=[y_3D,z_new];
        z_3D=[z_3D,yf];
        
%         hold on;
%         scatter3(x_3D,y_3D,z_3D);
%         axis equal;
%         hold off;
%            str=['alpha=',num2str(alpha/pi),'\pi. ','id=',num2str(id)];
%         title(str);
%        

    end
  
   
    
    
end
  [U,~]=alphavol([x_3D',y_3D'],2,0);
    dome_area=dome_area+U;
    


coord=[x_3D;y_3D;z_3D]';
coord=real(coord);
    
    hold on;
    [v,s]=alphavol(coord,1.5,1);
    shading interp;
    colormap copper;
    axis equal;
    axis tight;
    S=[S;s.bnd+size(PC,1)];
    PC=[PC;coord];




% xinv=diag([-1 1 1]);
% yinv=diag([1 -1 1]);
% S=[S;S+size(PC,1);S+2*size(PC,1);S+3*size(PC,1)];
Cratio=dome_area/(d0^2);
% PC=[PC; PC*xinv;PC*yinv;PC*xinv*yinv];


function diff=funfun(input,n,K,b_test,d0,H_test)

phi=input(1);
thetaa=input(2);
ksi=atan2(K*sin(pi/2+phi),cos(pi/2+phi))-pi/2; % vertical incidental angle
Theta=asin(sin(thetaa-ksi)/n)+ksi; % maximum vertical ray angle
thetac=asin(1/n); %critical angle 
theta0=pi-Theta-2*thetac; %exciting angle of parallel rays
d1=n*sin(theta0)*d0/sin(thetaa); %entrance aperture size by etuende conservation
b=d1/2*sqrt((tan(phi)^2+K)/(K*tan(phi)^2)); %semi minor axis;
diff(1)=b-b_test;
h=sign(phi)*b*(1-sqrt(K/(tan(phi)^2+K)));
Height=0.5*(d1+d0)*cot(Theta)+h;% trial height
diff(2)=Height-H_test;
end



