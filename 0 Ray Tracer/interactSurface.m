function [newray, powerout]=interactSurface(ray,surface,powerin)
%This function takes in a ray and a surface and outputs the ray that results
%from the interaction between the ray and the surface
powerout=powerin;
%calculate intersection
[intersection jout]=intersecSurface(surface, ray);
if length(jout)>1
    jout=jout(1);
    intersection=intersection(:,1);
end
if isempty(intersection)
    newray=createray(intersection,'dir',[0;0]);
    disp('no intersection')
    return
end

%determine if the ray is reflected or absorbed
switch surface.surftype
    case 'reflecting'
        %cosine of angle between ray and surface-parallel unit vectors
        phi=acos(dot(ray.pu,surface.tan)); 
        %decide between the cases: planeright vs. planeleft; v1, v2, v3, v4
        %The logic (1-truecase) is used to flip new_par when needed.
            new_par=sin(phi)*surface.norm*...
                (1-2*(surface.tan(1)>=0 && ray.pu(1)<0))*...
                (1-2*(surface.tan(1)<0 && ray.pu(1)<0))*...
                (1-2*(surface.tan(1)>0 && ray.pu(1)>0 &&...
                ray.pu'*surface.norm>0))*...
                (1-2*(surface.tan(1)<0 && ray.pu(1)>0 &&...
                ray.pu'*surface.norm>0));
            new_norm=cos(phi)*surface.tan;
        %new ray is comprised of two components relative to the surface:
        %parallel and normal; power is product of surface reflectivity and 
        %incident power
        newray=createray(intersection,'dir',new_par+new_norm,...
            surface.reflectivity*ray.power/100,ray.n);
        powerout.totalReflectionLoss=powerout.totalReflectionLoss+...
            abs(newray.power-ray.power);
    case 'absorbing'
        newray=createray(intersection,'dir',[0;0],0); %new ray has 0 power
        powerout.totalAbsorptionLoss=powerout.totalAbsorptionLoss+...
            abs(newray.power-ray.power,ray.n);
    case 'refracting'
        %angle between ray and surface-perpendicular unit vectors
        %Only work with 2d. Need to get normal vector at intersection
        normVecX=interp1((1:length(surface.tan)),surface.norm(:,1),jout);
        normVecY=interp1((1:length(surface.tan)),surface.norm(:,2),jout);
        if dot([normVecX;normVecY],ray.pu)>0
            surfaceNormal=[normVecX;normVecY];
        else
            surfaceNormal=[-normVecX;-normVecY];
        end
        %Define surface angle as absolute (0 radians would be || to [1;0])
        thetaSurfNormal=atan2(surfaceNormal(2),surfaceNormal(1));

        %Define angle of incoming ray
        thetaRay=atan2(ray.pu(2),ray.pu(1));

        %d_theta is the angle between the surface normal and the incoming ray
        d_theta=thetaSurfNormal-thetaRay; 
        
        %Use Snell's law to determine the new angle (relative to the
        %surface normal)
        n1=ray.n;
        [minDiff,I] = min([n1-surface.n1,n1-surface.n2]);
        if I==1
            n2=surface.n1;
        else
            n2=surface.n2;
        end
        newTheta = asin(n1*sin(d_theta)/n2);
        %Express angle relative to [1;0]
        newphi=thetaSurfNormal-newTheta;
        
        newray=createray(intersection,'dir',[cos(newphi);sin(newphi)],...
            ray.power,n2);
    otherwise
        beep
        disp('invalid surface type')
        return
end
powerout.Remaining=powerout.totalIncident-powerout.totalReflectionLoss-...
    powerout.totalAbsorptionLoss;
        