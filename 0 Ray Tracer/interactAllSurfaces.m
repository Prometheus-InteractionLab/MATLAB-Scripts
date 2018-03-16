function [newray, powerout,surfaceNum]=...
    interactAllSurfaces(ray,surfaces,powerin)
%This function takes in a ray and a vector of surfaces and outputs the ray 
%that results from the interaction between the ray and the surface. The ray
%can be reflected, absorbed, or refracted, depending on the surface type.
powerout=powerin;
%calculate intersection (yes, this has been done before this function is
%called, but the program is modular and this is a bit easier to work with).
[surfaceNum,intersectionOut,jout]=intersecAllSurfaces(surfaces, ray);
%If the ray doesn't intersect with any surfaces, "kill" it by giving it a
%direction that is the zero vector
if isempty(intersectionOut)
    newray=createray(intersectionOut,'dir',[0;0]);
    disp('no intersection')
    return
else
    intersection=intersectionOut; %point of intersection
    surface=surfaces{surfaceNum}; %surface with which the ray interacts
end

%determine if the ray is reflected, refracted, or absorbed
switch surface.surftype
    case 'reflecting'
        %cosine of angle between ray and surface-parallel unit vectors
        phi=acos(dot(ray.pu,surface.tan)); 
        %Compute the new angle - this is a bit complicated for the general
        %case but the logic here works.
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
        %the ray is simply absorbed. Replace it with a ray starting at the
        %intersection and having a direction that is the zero vector and a
        %power of 0. This is used in simulations to "kill" rays so that the
        %surviving rays can more easily be identified.
        newray=createray(intersection,'dir',[0;0],0,1); %new ray has 0 power
%         powerout.totalAbsorptionLoss=powerout.totalAbsorptionLoss+...
%             abs(newray.power-ray.power,ray.n);
    case 'refracting'
        %We need the angle between the ray and surface-perpendicular unit
        %vectors. This currently works in 2D.
        %Get the surface normal vector at the intersection
        normVecX=interp1((1:length(surface.tan)),surface.norm(:,1),jout);
        normVecY=interp1((1:length(surface.tan)),surface.norm(:,2),jout);
        if dot([normVecX;normVecY],ray.pu)>0
            surfaceNormal=[normVecX;normVecY];
        else
            surfaceNormal=[-normVecX;-normVecY];
        end
        %Define surface angle as absolute (0 radians would be || to [1;0])
        thetaSurfNormal=atan2(surfaceNormal(2),surfaceNormal(1));
        %Define the angle of the incoming ray
        thetaRay=atan2(ray.pu(2),ray.pu(1));
        %d_theta is the angle between the surface normal and the incoming ray
        d_theta=thetaSurfNormal-thetaRay; 
        %Use Snell's law to determine the new angle (relative to the
        %surface normal)
        n1=ray.n;
        %n2 is given by the surface. It isn't known which side the ray is
        %on, but the new index of refraction should be as different as
        %possible from the current one. If there's no difference in the
        %index of refraction on either side of the surface, then it doesn't
        %matter which one is the new n
        [maxDiff,I] = max([abs(n1-surface.n1),abs(n1-surface.n2)]);
        if I==1
            n2=surface.n1;
        else
            n2=surface.n2;
        end
        newTheta = asin(n1*sin(d_theta)/n2);
        if imag(newTheta)~=0;
            %internal reflection - the index of refraction of the ray stays
            %the same and the ray is reflected within its current object
            newphi=thetaSurfNormal+pi+d_theta;
            newray=createray(intersection,'dir',[cos(newphi);sin(newphi)],...
                ray.power,ray.n);
        else
            %The ray is refracted and takes on a different direction and
            %index of refraction
            %Express angle relative to [1;0]
            newphi=thetaSurfNormal-newTheta;
            %newphi*180/pi
            newray=createray(intersection,'dir',[cos(newphi);sin(newphi)],...
                ray.power,n2);
        end

        
    otherwise
        %This will throw an error and stop execution since the function
        %doesn't define newray in this case.
        beep
        fprintf('invalid surface type: %s\n',surface.surftype);
        return
end
powerout.Remaining=powerout.totalIncident-powerout.totalReflectionLoss-...
    powerout.totalAbsorptionLoss;
        