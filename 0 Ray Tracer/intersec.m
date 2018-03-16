function intersection=intersec(plane, ray)
%for a given plane and ray, calculate the intersection if there is one

P1=plane.center-plane.radius*plane.parallel;    %startpoint of plane
P2=plane.center+plane.radius*plane.parallel;
if P1(2)>P2(2)  %make sure parallel to plane points up
   temp=P1;
   P1=P2;
   P2=temp;
end

S=ray.t_ray(0);                                 %startpoint of ray

t_intersec=dot(P1-S,plane.norm)/dot(ray.pu,plane.norm);
temp=ray.t_ray(t_intersec);
if t_intersec<0 %in forward direction of ray?
    intersection = [];
%     disp('plane behind ray')
elseif (temp(2)<P1(2) || temp(2)>P2(2)) ||...
        (temp(1) < min(P1(1),P2(1)) || temp(1) > max(P1(1),P2(1)))
    %intersection out of bounds?
    intersection = [];
%     disp('ray doesn''t intersect plane')
else
    intersection=temp;
end


