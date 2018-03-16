function surface=createSurface(x,y,order,dt,surftype,n1,n2)
    %This function creates a 2D surface object. The inputs are
    %   x,y: coordinates of points that make up the surface
    %   order: what is the order of the spline that will be fitted to x,y
    %   dt: how fine should the resolution of the spline be? smaller
    %       numbers here mean better resolution but slower computation
    %   surftype: does the surface reflect, refract, or absorb? refracting
    %             surfaces can exhibit total internal reflection with the
    %             non-sequential ray tracer
    %   n1: the index of refraction on one side
    %   n2: the index of refraction on the other side
    % Note that a ray is defined as having an index of refraction. When it
    % passes through a surface and is refracted, it is assumed that the ray
    % has an index of refraction equal to either n1 or n2. The new index of
    % refraction is then assumed to be n2 or n1.
        
    %set parameter t and the resampled ts
    t=linspace(0,1,length(x));
    ts=0:dt:1;
    disp('Makin'' splines for the surface...')
    val = deboor(t,[x,y],ts,order);
    xs=val(:,1);
    ys=val(:,2);
%     xs = (x(1):dt:x(end))';
%     ys = spline(x,y,xs);
    
    
    %Get tangent vectors (parallel to surface). This could easily be
    %extended if a z coordinate were specified
    tanvec=[gradient(xs),gradient(ys)]; 
%     figure, plot(tanvec(:,1),tanvec(:,2))
    %unit tangent vector
    U_tanvec(:,1)=tanvec(:,1)./sqrt(tanvec(:,1).^2+tanvec(:,2).^2);
    U_tanvec(:,2)=tanvec(:,2)./sqrt(tanvec(:,1).^2+tanvec(:,2).^2);
    U_tanvec(:,3)=zeros(size(U_tanvec,1),1);

    %Get normal vector - cross product between tangential vector and k-hat
    %vector. This works very well in 2D, but would have to be changed in 3D
    U_normvec=zeros(size(U_tanvec));
    for r=1:size(U_normvec,1);
        U_normvec(r,:)=cross(U_tanvec(r,:),[0, 0, -1]);
    end
    
    %Now create the surface object
    surface.x=xs;
    surface.y=ys;
    surface.tan=U_tanvec;
    surface.norm=U_normvec;
    surface.surftype=surftype;
    surface.reflectivity=0;
    surface.n1=n1; %What was the index of refraction?
    surface.n2=n2; %What is the index of refraction on the other side?
    
%     figure
%     plot(x,y,'.')
%     hold on
%     plot(xs,ys,'r','LineWidth',1)
%     axis equal
%     quiver(xs,ys,U_normvec(:,1),U_normvec(:,2))

function val = deboor(T,p,y,order)
% function val = DEBOOR(T,p,y,order)
% val = DEBOOR(t,[x,y],ts,4);
% INPUT:  T     Stützstellen
%         p     Kontrollpunkte (nx2-Matrix)
%         y     Auswertungspunkte (Spaltenvektor)
%         order Spline-Ordnung
%
% OUTPUT: val   Werte des B-Splines an y (mx2-Matrix)
%
% Date:   2007-11-27
% Author: Jonas Ballani

m = size(p,1); %how many control points?
n = length(y); %how many points are to be output?
X = zeros(order,order);
Y = zeros(order,order);
a = T(1);
b = T(end);
T = [ones(1,order-1)*a,T,ones(1,order-1)*b];

%loop through each of the points that are to be output
for l = 1:n
    t0 = y(l); %current value of parameter t
    id = find(t0 >= T);
    k = id(end);
		if (k > m)
			return;
		end
    X(:,1) = p(k-order+1:k,1);
    Y(:,1) = p(k-order+1:k,2);

    for i = 2:order
        for j = i:order
            num = t0-T(k-order+j);
            if num == 0
                weight = 0;
            else
                s = T(k+j-i+1)-T(k-order+j);
                weight = num/s;
            end
            X(j,i) = (1-weight)*X(j-1,i-1) + weight*X(j,i-1);
            Y(j,i) = (1-weight)*Y(j-1,i-1) + weight*Y(j,i-1);
        end
    end
    val(l,1) = X(order,order);
    val(l,2) = Y(order,order);
end