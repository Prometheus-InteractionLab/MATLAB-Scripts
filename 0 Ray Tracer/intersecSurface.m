function [intersection,jout]=intersecSurface(surface, ray)
%for a given surface and ray, calculate the intersection if there is one


P0=ray.t_ray(0);                                 %startpoint of ray
P1=ray.t_ray(1000);                              %second point on ray
%1 is the ray and 2 is the surface

[x0,y0,iout,jout] = computeIntersections([P0(1),P1(1)],[P0(2),P1(2)],...
    surface.x,surface.y,true);
if length(jout)>1
    intersection=zeros(2,length(jout));
    for idx=1:length(jout)
       intersection(:,idx)=[x0(idx);y0(idx)]; 
    end
else
    intersection = [x0;y0];
end
% t_intersec=dot(P1-S,plane.norm)/dot(ray.pu,plane.norm);
% t_intersec=dot(P1-S,plane.norm)/dot(ray.pu,plane.norm);
% temp=ray.t_ray(t_intersec);
% if t_intersec<0 %in forward direction of ray?
%     intersection = [];
% %     disp('plane behind ray')
% elseif (temp(2)<P1(2) || temp(2)>P2(2)) ||...
%         (temp(1) < min(P1(1),P2(1)) || temp(1) > max(P1(1),P2(1)))
%     %intersection out of bounds?
%     intersection = [];
% %     disp('ray doesn''t intersect plane')
% else
%     intersection=temp;
% end


