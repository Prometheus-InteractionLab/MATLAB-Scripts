function [newray, powerout]=interact(ray,plane,powerin)
%This function takes in a ray and a plane and outputs the ray that results
%from the interaction between the ray and the plane
powerout=powerin;
%calculate intersection
intersection=intersec(plane, ray);
if isempty(intersection)
    newray=createray(intersection,'dir',[0;0]);
    disp('no intersection')
    return
end

%determine if the ray is reflected or absorbed
switch plane.surftype
    case 'reflecting'
        %cosine of angle between ray and plane-parallel unit vectors
        phi=acos(dot(ray.pu,plane.parallel)); 
        %decide between the cases: planeright vs. planeleft; v1, v2, v3, v4
        %The logic (1-truecase) is used to flip new_par when needed.
            new_par=sin(phi)*plane.norm*...
                (1-2*(plane.parallel(1)>=0 && ray.pu(1)<0))*...
                (1-2*(plane.parallel(1)<0 && ray.pu(1)<0))*...
                (1-2*(plane.parallel(1)>0 && ray.pu(1)>0 &&...
                ray.pu'*plane.norm>0))*...
                (1-2*(plane.parallel(1)<0 && ray.pu(1)>0 &&...
                ray.pu'*plane.norm>0));
            new_norm=cos(phi)*plane.parallel;
        %new ray is comprised of two components relative to the plane:
        %parallel and normal; power is product of plane reflectivity and 
        %incident power
        newray=createray(intersection,'dir',new_par+new_norm,...
            plane.reflectivity*ray.power/100);
        powerout.totalReflectionLoss=powerout.totalReflectionLoss+...
            abs(newray.power-ray.power);
    case 'absorbing'
        newray=createray(intersection,'dir',[0;0],0); %new ray has 0 power
        powerout.totalAbsorptionLoss=powerout.totalAbsorptionLoss+...
            abs(newray.power-ray.power);
    case 'refracting'
        %angle between ray and plane-perpendicular unit vectors
        phi=-acos(dot(ray.pu,plane.norm)); 
        newphi = asin(plane.nOld*sin(phi)/plane.nNew);
        newray=createray(intersection,'dir',[cos(newphi);sin(newphi)],...
            ray.power);
    otherwise
        beep
        disp('invalid surface type')
        return
end
powerout.Remaining=powerout.totalIncident-powerout.totalReflectionLoss-...
    powerout.totalAbsorptionLoss;
        