n=1.41; %refractive index
phi=pi/12; %predefined arc angle at x-axis
thetaa1=21/180*pi; % acceptance angle at x-axis
thetaa2=pi/6; % acceptance angle at y-axis
d00=5; %size of the bottom square aperture
Theta=asin(sin(thetaa1-phi)/n)+phi; % maximum vertical ray angle
thetac=asin(1/n); %critical angle 
theta0=pi-Theta-2*thetac; %exciting angle of parallel rays
d1=n*sin(theta0)*d00/sin(thetaa1); %entrance aperture size by etuende conservation
R=d1/2/sin(phi); % radius
h=R*(1-cos(phi)); %height of the dome
H0=0.5*(d1+d00)*cot(Theta)+h; %height of the collector.
PC=[];
S=[];
dome_area=0;
x_top=[];
y_top=[];
z_top=[];
x_bracket=[];
z_bracket=[];
x3d=[];
y3d=[];
z3d=[];
for range=0:1
  
    for alpha=(range)*pi/4:pi/4/16:(range+1)*pi/4 %rotation angle
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
        d0=d00*sqrt(1+tan(alpha)^2); % exit aperture
        if(alpha>pi/4)
            d0=d00*sqrt(1+cot(alpha)^2);
%             thetaa=thetaa2;
        end
        
        fun=@(x) f(x,n,R,d0,H0); %solve for the arc angle
        sols=fsolve(fun,[phi,thetaa1]);
        phi=sols(1);
        thetaa=sols(2);
        
        theta=linspace(-phi,phi,160)';
        thetap=asin(sin(theta+thetaa)/n)-theta;
        Theta=asin(sin(thetaa-phi)/n)+phi; % maximum vertical ray angle
        thetac=asin(1/n);
        theta0=pi-Theta-2*thetac;
        
        d1=n*sin(theta0)*d0/sin(thetaa);
        R=d1/2/sin(phi); % trial radius
        l1=2*R*sin((theta+thetaa)/2).^2; %segment 1 of the optical path
        h=R*(1-cos(phi));
        H0=0.5*(d1+d0)*cot(Theta)+h;
        
        C=2*R*sin((thetaa+phi)/2)^2+n*(d1+d0)/(2*sin(Theta)); % constant optical path
        h=R*(1-cos(phi));
        
        
        y0=R*(1-cos(theta));
        a=(C-l1)/n;
        b=R*sin(theta)+d0/2;
        c=H0-y0;
        l2=(a.^2-b.^2-c.^2)./(2*(a+b.*sin(thetap)-c.*cos(thetap))); %segment 2 for rays hitting the same point
        x=R*sin(theta)+l2.*sin(thetap);
        z=H0-l2.*cos(thetap)-y0;
              
        
        l2=(a-b*sin(theta0)-c*cos(theta0))./(2*sin((theta0+thetap)/2).^2); %segment 2 for rays exiting parallelly 
        xp=R*sin(theta)+l2.*sin(thetap);
        zp=H0-l2.*cos(thetap)-y0;
        
        d0_new=2*xp(1); %exiting aperture determined by the side profile, should be equal to d0
               
        [X,XP]=meshgrid(x,xp); %findin the intersection of 2 choices of segment 2
        [Z,ZP]=meshgrid(z,zp);
        D=(X-XP).^2+(Z-ZP).^2;
        Dmin=min(min(D));
        [~,t]=find(Dmin==D);
        z0=z(t); % the intersection 
        x_dome=linspace(d1/2/16,d1/2,16)';
        z_dome=sign(R)*sqrt(R^2-x_dome.^2)+H0-R;
        x_top=[x_top,x_dome*cos(alpha)];
        y_top=[y_top,x_dome*sin(alpha)];
        z_top=[z_top,z_dome];
        xf=(z>z0).*x+(z<=z0).*xp; %The side profile at a specific angle, including the dome
        zf=(z>z0).*z+(z<=z0).*zp;
        
            
        
        if alpha<eps||(abs(alpha-pi/4)<eps)||(abs(alpha-pi/2)<eps)
            %generate the bracket plane at 0, pi/4
            for k=1:length(zf)
                if zf(k)<H0
                    x_temp=linspace(0,d1,20);
                    x_side=x_temp(x_temp<=xf(k));
                    x_bracket=[x_bracket,x_side];
                    z_bracket=[z_bracket,repmat(zf(k),1,length(x_side))];
                end
            end
        end
        
        
        x_new=xf*cos(alpha); %Coordinates Transformation
        y_new=xf*sin(alpha);
        
        x3d=[x3d,x_new]; %Append the entire point sets with the side profile
        y3d=[y3d,y_new];
        z3d=[z3d,zf];
       

    end
  
    
    

%     [v,s]=alphavol(coordinates,8,1);
%   
%     shading interp;
%     colormap copper;
%     axis equal;
%     axis tight;
%     S=[S;s.bnd+size(PC,1)];
%     PC=[PC;coordinates];
    
    
end
[a,~]=alphavol([reshape(x_top,[],1),reshape(y_top,[],1)],3,1);
dome_area=dome_area+a;


xinv=diag([-1 1 1]);
yinv=diag([1 -1 1]);
% S=[S;S+size(PC,1);S+2*size(PC,1);S+3*size(PC,1)];
dome_area=dome_area*4;
C_ratio=dome_area/(d0^2);
% PC=[PC; PC*xinv;PC*yinv;PC*xinv*yinv];

figure;
x3d_all=[x3d,-flip(x3d,2),-x3d,flip(x3d,2)];
y3d_all=[y3d,flip(y3d,2),-y3d,-flip(y3d,2)];
z3d_all=[z3d,z3d,z3d,z3d];

figure;surf(x3d_all,y3d_all,z3d_all);shading interp;
colormap(copper);

axis equal;



% % V1=[PC(S(:,1),1),PC(S(:,1),2),PC(S(:,1),3)];
% % V2=[PC(S(:,2),1),PC(S(:,2),2),PC(S(:,2),3)];
% % V3=[PC(S(:,3),1),PC(S(:,3),2),PC(S(:,3),3)];
% % l1=V2-V1;
% % l2=V3-V2;
% faceNorm=cross(l1,l2,2);
% idx=find(faceNorm(:,3)==0);
% S(idx,:)=[];
% FnV=struct('vertices',PC,'faces',S);
% stlwrite('DTIRC0211_pi8.STL',FnV);

function diff=f(input,n,R_trial,d0,H_trial)

phi=input(1);
thetaa=input(2);

Theta=asin(sin(thetaa-phi)/n)+phi;

thetac=asin(1/n); %critical angle 
theta0=pi-Theta-2*thetac; %exciting angle of parallel rays
d1=n*sin(theta0)*d0/sin(thetaa); %entrance aperture size by etuende conservation
R=d1/2/sin(phi); % radius
diff(1)=R-R_trial;

h=R*(1-cos(phi));

Height=0.5*(d1+d0)*cot(Theta)+h;% trial height
diff(2)=Height-H_trial;
end



