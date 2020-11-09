PC=[]; %The coordinates of the pointCloud to be generated. 
dome_area=0; % The xy plane projected area of the top surface
x_top=[]; % top surface coordinates on the computation plane.
y_top=[];
n=1.41;
thetaa1=pi/7.2; % acceptance angle
pAngle=-pi/15; % ramp angle
d00=5.0; %
d0=d00;
d1=15; % attepted d1
thetaa=thetaa1; 
Theta=asin(sin(thetaa1-pAngle)/n)+pAngle; %Vertical inclination of the 
%extreme ray hitting the left corner of the top surface 
H0=0.5*(d1+d0)*cot(Theta)+d1/2*tan(pAngle); %The height at the center, 
%consistent for side profiles of all angles to form a continuous top dome shape

%Solve for H and d1 from all the input parameters
func=@(x) DTIRCtrialParameters(n,d0,x(1),pAngle,thetaa1,x(2));
sol=fsolve(func,[d1,H0]);
d1=sol(1);
H0=sol(2);
H=0.5*(d1+d0)*cot(Theta);

%generate the x,y coordinates
[diff,points,idx]=DTIRCtrialParameters(n,d0,sol(1),pAngle,thetaa1,sol(2));
x=points{1}(:,1);
z=points{1}(:,2);
xp=points{2}(:,1);
zp=points{2}(:,2);
[X,XP]=meshgrid(x,xp); %find the intersection of 2 choices of segment 2
[Z,ZP]=meshgrid(z,zp);
D=(X-XP).^2+(Z-ZP).^2;
Dmin=min(min(D));
[~,t]=find(Dmin==D);
z0=z(t); % the intersection
%generate xf, zf final coordinates
xf=[(z>z0).*x+(z<=z0).*xp];
zf=[(z>z0).*z+(z<=z0).*zp];
% xfi=xf;
% zfi=zf;

%interpolation of the sidewall
zfi_upper=linspace(zf(idx),zf(end),260-idx+1);
xfi_upper=interp1(zf(idx:end),xf(idx:end),zfi_upper,'spline','extrap');
xfi=[xf(1:idx-1);xfi_upper'];
zfi=[zf(1:idx-1);zfi_upper'];

%side profile coordinates stored in surface matrix
x3d=xfi;
z3d=zfi;
y3d=zeros(1,length(x3d))';

%generate the top coordinates
xtop=linspace(0.02,xf(end),20)';
ytop=zeros(size(xtop,1),size(xtop,2));
xfztop=H0-sqrt(xtop.^2+ytop.^2)*tan(pAngle);

%Generate the side profiles for all azimuthal angles
for alpha=pi/4/16:pi/4/16:pi/2 %rotation angle
    d0=d00*sqrt(1+tan(alpha)^2); % exit aperture
    if(alpha>pi/4)
        d0=d00*sqrt(1+cot(alpha)^2);
    end
    func=@(x) DTIRCtrialParameters(n,d0,x(1),pAngle,x(2),H0);
    sol=fsolve(func,[d1,thetaa]);
    [diff,points,idx]=DTIRCtrialParameters(n,d0,sol(1),pAngle,sol(2),H0);
    fprintf('alpha=%f, d0=%f, \n.', alpha/pi*180,d0);
    d1=sol(1);
    thetaa=sol(2)；
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
%     xfi=xf;
%     zfi=zf;
    zfi_upper=linspace(zf(idx),zf(end),260-idx+1);
    xfi_upper=interp1(zf(idx:end),xf(idx:end),zfi_upper,'spline','extrap');
    xfi=[xf(1:idx-1);xfi_upper'];
    zfi=[zf(1:idx-1);zfi_upper'];
    x3d=[x3d,xfi*cos(alpha)];
    z3d=[z3d,zfi];
    y3d=[y3d,xfi*sin(alpha)];
    r_top=linspace(0.02,xf(end),20)';
    xtop_new=r_top*cos(alpha);
    ytop_new=r_top*sin(alpha);
    ztop_new=H0-sqrt(xtop_new.^2+ytop_new.^2)*tan(pAngle);
    xtop=[xtop,xtop_new];
    ytop=[ytop,ytop_new];
    ztop=[ztop,ztop_new];
end
%Organize the point cloud
x3d=[fliplr(x3d),x3d,fliplr(-x3d),-x3d];
y3d=[-fliplr(y3d),y3d,fliplr(y3d),-y3d];
z3d=[fliplr(z3d),z3d,fliplr(z3d),z3d];



%calculate concentration ratio
[S,~]=alphavol([reshape(xtop,[],1),reshape(ytop,[],1)],1,1);
C_ratio=4*S/d00^2;
%process the data for stl file conversion
P1=1:length(x3d)/2;
P2=(length(x3d)/2+1):length(x3d);
[F,V]=surf2patch(x3d,y3d,z3d,'triangles');
[F1,V1]=surf2patch(x3d(P1,:),y3d(P1,:),z3d(P1,:),'triangles');
[F2,V2]=surf2patch(x3d(P2,:),y3d(P2,:),z3d(P2,:),'triangles');

%The function to solve for a consistent set of parameters
function [diff,points,idx]=DTIRCtrialParameters(n,d0,d1,pAngle,thetaa,H0_trial)
Theta=asin(sin(thetaa-pAngle)/n)+pAngle; % maximum vertical ray angle
thetac=asin(1/n);  %critical angle of the material
theta0=pi-Theta-2*thetac; 
%Angle of reflected ray from the extreme ray hitting the left corner of the top surface

xd=linspace(-d1/2,d1/2,160); %the locations of rays hitting the top plane
%the REAL locations of rays hitting the top surfaces.
xdpl=(xd+d1/2)/(1-tan(-pAngle)*tan(thetaa))-d1/2;
ydpl=-(xd+d1/2)*tan(-pAngle)/(1-tan(-pAngle)*tan(thetaa));
xdpr=d1/2-(d1/2-xd)/(1+tan(-pAngle)*tan(thetaa));
ydpr=-(d1/2-xd)*tan(-pAngle)/(1+tan(-pAngle)*tan(thetaa));
xdc=d1/2*tan(pAngle)*tan(thetaa);
%the location of the ray (in terms of xd) hitting exactly the center of the top surface
idx=sum(xd<xdc);
xdp=[xdpl(1:idx),xdpr(idx+1:end)];
ydp=[ydpl(1:idx),ydpr(idx+1:end)];

pA=sign(xdp)*(-pAngle);
thetap=asin(sin(thetaa-pA)/n)+pA;%Vertical inclination angles of all rays
l1=((xdp+d1/2)*tan(thetaa)-ydp)*cos(thetaa);
H=0.5*(d1+d0)*cot(Theta);
C=d1*sin(thetaa)+n*(d1+d0)/(2*sin(Theta)); 
%The total OPL: constant optical path for ray hitting the left corner (coming from left)
dC=-0.6; %OPL difference across the center of the top surface
deltaC=[dC*ones(1,idx),zeros(1,length(xd)-idx)];
%deltaC=flip(cumtrapz(d1/2-flip(xd),(n*flip(sin(thetap))-sin(thetaa))));
C=deltaC+C;
a=(C-l1)/n;
b=xdp+d0/2;
c=H+ydp;
l2=(a.^2-b.^2-c.^2)./(2*(a+b.*sin(thetap)-c.*cos(thetap))); %segment 2 for rays hitting the same point
x=xdp+l2.*sin(thetap);%Side profile coordinates generated from knowing l2;
y=H+ydp-l2.*cos(thetap);
l2=(a-b*sin(theta0)-c*cos(theta0))./(2*sin((theta0+thetap)/2).^2); %segment 2 for rays exiting parallelly
xp=xdp+l2.*sin(thetap);
yp=H+ydp-l2.*cos(thetap);
H0=H+d1/2*tan(pAngle);
d0_new=2*xp(1); %exiting aperture determined by the side profile, should be equal to d0
diff(1)=d0-d0_new;
diff(2)=H0-H0_trial;
points={[x',y'],[xp',yp']};
end
