function [surfaceNum,intersectionOut,jout]=...
    intersecAllSurfaces(surfaces, ray)
%given the set of surfaces and a single ray, return the closest
%intersection and the index of the surface

P0=ray.t_ray(0);                    %startpoint of ray
P1=ray.t_ray(10000);                %second point on ray (10 m)

rMin=inf;  %initialize the minimum radius as infinite
intersectionOut=[]; %Initialize as no intersection
jout = []; %Initialize as no intersection
surfaceNum=[];
%loop through all of the surfaces and find the intersection closest to P0
for s=1:length(surfaces)
    surface=surfaces{s};
    %The routine to compute the intersection of a line and a surface was
    %obtained from the MATLAB file exchange. It seems to be robust and
    %effective, but not super fast. This would need to be changed for a 3D
    %surface
    [x0,y0,iout,j] = computeIntersections([P0(1),P1(1)],[P0(2),P1(2)],...
        surface.x,surface.y,true); %#ok<ASGLU>
    if length(j)>1
        %more than one intersection
        r=zeros(length(j),1);
        for idx=1:length(j)
           r(idx)=sqrt((P0(1)-x0(idx))^2+(P0(2)-y0(idx))^2);
           if r(idx)<1E-8 %it is the start point
               r(idx)=inf;
           end
        end
        %For the given surface, determine the minimum radius
        [smallestRforSurface, I] = min(r);
        intersection = [x0(I);y0(I)];
        j=j(I);
    else
        
        intersection = [x0;y0];
        
        if isempty(intersection)
            %Intersection never occurs
            smallestRforSurface=inf;
        else
            r=sqrt((P0(1)-x0)^2+(P0(2)-y0)^2);
            %For the given surface, determine the minimum radius
            if r < 1E-8
                %Assume that this is the current surface - ignore it
                smallestRforSurface=inf;
            else
                smallestRforSurface = r;
            end
        end
    end
    %see if smallest radius for this surface is smaller than the current
    %smallest radius
    if smallestRforSurface < rMin
       %the surface number needs to be updated, as does intersection and
       %jout
       rMin = smallestRforSurface;
       surfaceNum = s;
       intersectionOut = intersection;
       jout = j;
    end
end



