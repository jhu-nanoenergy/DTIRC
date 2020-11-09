PC=[];
dome_area=0;
x_top=[];
y_top=[];
n=1.41;
thetaa1=pi/8;
pAngle=0;
d00=5.0;
d0=d00;
d1=15;
thetaa=thetaa1;
Theta=asin(sin(thetaa1-pAngle)/n)+pAngle;
H=0.5*(d1+d0)*cot(Theta);


func=@(x) DTIRCtrialParameters(n,d0,x(1),pAngle,thetaa1,x(2));
sol=fsolve(func,[d1,H]);
[diff,points]=DTIRCtrialParameters(n,d0,sol(1),pAngle,thetaa1,sol(2));
d1=sol(1);
H=sol(2);


x=points{1}(:,1);
z=points{1}(:,2);
xp=points{2}(:,1);
zp=points{2}(:,2);
[X,XP]=meshgrid(x,xp); %findin the intersection of 2 choices of segment 2
[Z,ZP]=meshgrid(z,zp);
D=(X-XP).^2+(Z-ZP).^2;
Dmin=min(min(D));
[~,t]=find(Dmin==D);
z0=z(t); % the intersection

xf=[(z>z0).*x+(z<=z0).*xp];
zf=[(z>z0).*z+(z<=z0).*zp];
% xx=[xp(1:t);x(51:100)];
% zz=[zp(1:t);z(51:100)];
% zf=linspace(0,H,200)';
% xf=interp1(zz,xx,zf,'spline');

x3d=xf;
z3d=zf;
y3d=zeros(1,length(x3d))';




for alpha=pi/4/16:pi/4/16:pi/4 %rotation angle
    d0=d00*sqrt(1+tan(alpha)^2); % exit aperture
    if(alpha>pi/4)
        d0=d00*sqrt(1+cot(alpha)^2);
        %             thetaa=thetaa2;
    end
    
    
    func=@(x) DTIRCtrialParameters(n,d0,x(1),pAngle,x(2),H);
    sol=fsolve(func,[d1,thetaa]);
    [diff,points]=DTIRCtrialParameters(n,d0,sol(1),pAngle,sol(2),H);
    fprintf('alpha=%f, d0=%f, \n.', alpha/pi*180,d0);
    
    
    d1=sol(1);
    thetaa=sol(2);
    
    
    x=points{1}(:,1);
    z=points{1}(:,2);
    xp=points{2}(:,1);
    zp=points{2}(:,2);
    [X,XP]=meshgrid(x,xp); %findin the intersection of 2 choices of segment 2
    [Z,ZP]=meshgrid(z,zp);
    D=(X-XP).^2+(Z-ZP).^2;
    Dmin=min(min(D));
    [~,t]=find(Dmin==D,1);
    z0=z(t); % the intersection
    
    xf=(z>z0).*x+(z<=z0).*xp;
    zf=(z>z0).*z+(z<=z0).*zp;
    
    
%     xx=[xp(1:t);x(51:100)];
% zz=[zp(1:t);z(51:100)];
% zf=linspace(0,H,200)';
% xf=interp1(zz,xx,zf,'spline');
    
    x3d=[x3d,xf*cos(alpha)];
    z3d=[z3d,zf];
    y3d=[y3d,xf*sin(alpha)];
    
end

[Xtop,Ytop]=meshgrid(linspace(0,d1/2,100));
xtop=reshape(Xtop,1,numel(Xtop));
ytop=reshape(Ytop,1,numel(Ytop));

for k=1:length(xtop)
    Ds=(xtop(k)-x3d(end,:)).^2+(ytop(k)-y3d(end,:)).^2;
    [~,ii]=find(Ds==min(Ds),1);
    if xtop(k)^2+ytop(k)^2>=0.98*(x3d(end,ii)^2+y3d(end,ii)^2)
        xtop(k)=NaN;
        ytop(k)=NaN;
    end
end

xtop=xtop(~isnan(xtop));
ytop=ytop(~isnan(ytop));
ztop=H*ones(1,length(xtop));

[aaa,bbb]=alphavol([[xtop,x3d(end,:)];[ytop,y3d(end,:)]]',1,1);
C_ratio=4*aaa/d00^2;


x3d=[fliplr(x3d),x3d];
y3d=[-fliplr(y3d),y3d];
z3d=[fliplr(z3d),z3d];
[F,V]=surf2patch(x3d,y3d,z3d,'triangles');
P1=1:length(x3d)/2;
P2=(length(x3d)/2+1):length(x3d);
[F1,V1]=surf2patch(x3d(P1,:),y3d(P1,:),z3d(P1,:),'triangles');
[F2,V2]=surf2patch(x3d(P2,:),y3d(P2,:),z3d(P2,:),'triangles');


function [diff,points]=DTIRCtrialParameters(n,d0,d1,pAngle,thetaa,H_trial)
Theta=asin(sin(thetaa-pAngle)/n)+pAngle; % maximum vertical ray angle
thetac=asin(1/n);
theta0=pi-Theta-2*thetac;
xd=linspace(-d1/2,d1/2,200);
pA=sign(xd).*abs(2*xd/d1).^0.001*(-pAngle);
thetap=asin(sin(thetaa-pA)/n)+pA;
l1=(xd+d1/2)*sin(thetaa); %segment 1 of the optical path
H=0.5*(d1+d0)*cot(Theta);
C=d1*sin(thetaa)+n*(d1+d0)/(2*sin(Theta));
%constant optical path for ray hitting the left corner (coming from left)
discLoc=length(xd)/2+1;
deltaC_right=flip(cumtrapz(d1/2-flip(xd(discLoc:end)),(n*flip(sin(thetap(discLoc:end)))-sin(thetaa))));
deltaC_left=+flip(cumtrapz(d1/2-flip(xd(1:discLoc-1)),(n*flip(sin(thetap(1:discLoc-1)))-sin(thetaa))));
%Difference in path length caused by Fresnel surface



deltaC=[deltaC_left,deltaC_right];%
deltaC=flip(cumtrapz(d1/2-flip(xd),(n*flip(sin(thetap))-sin(thetaa))));
C=deltaC+C;
a=(C-l1)/n;
b=xd+d0/2;
c=H;
l2=(a.^2-b.^2-c.^2)./(2*(a+b.*sin(thetap)-c.*cos(thetap))); %segment 2 for rays hitting the same point
x=xd+l2.*sin(thetap);
y=H-l2.*cos(thetap);


l2=(a-b*sin(theta0)-c*cos(theta0))./(2*sin((theta0+thetap)/2).^2); %segment 2 for rays exiting parallelly
xp=xd+l2.*sin(thetap);
yp=H-l2.*cos(thetap);

d0_new=2*xp(1); %exiting aperture determined by the side profile, should be equal to d0
diff(1)=d0-d0_new;
diff(2)=H-H_trial;
points={[x',y'],[xp',yp']};
end
