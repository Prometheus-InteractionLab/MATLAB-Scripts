function [power intersections rays prevrays surfaceNums]=...
    nonSequential_Surfacetrace(surfaces,rays,showplot)
%non-sequential ray tracer with reflection, refraction, and absorptions

%Initialize intersections cell - this is a cell of 2 x n matrices, where n
%is the number of intersections for a given ray. The first point is
%the starting point of each ray, and the following points are added to each
%row of the intersections cell as the simulation progresses. Each row of
%the cell corresponds to a ray.
intersections=cell(length(rays),1);
%Define total incident power at the beginning. At the moment this isn't
%really used, rather it's kept around for compatibility for future versions
power.totalIncident=0;
power.totalReflectionLoss=0;
power.totalAbsorptionLoss=0;

for r=1:length(rays)
    power.totalIncident=power.totalIncident+rays{r}.power;
end
power.Remaining=power.totalIncident-power.totalReflectionLoss-...
    power.totalAbsorptionLoss;


%set up initial points
for r=1:length(rays)
    intersections{r}=rays{r}.t_ray(0);
end
%plot all surfaces
if showplot
    figure;
    for s=1:length(surfaces)
       plot(surfaces{s}.x,surfaces{s}.y,'k-','Linewidth',2); 
       hold all
    end
end

prevrays=cell(size(rays));
surfaceNums=zeros(length(rays),1);
%Loop through each ray and trace it until its children don't intersect with
%any more surfaces
for r=1:length(rays)
    if mod(r,100)==0 || r==1
        fprintf('Calculation: Ray %2.0f/%2.0f\n',r,length(rays));
    end
    pause(0.01)
    ray=rays{r}; 
    %look at all surfaces this ray interacts with. Take the
    %closest intersection to start
    [surfaceNum,intersectionOut,jout]=intersecAllSurfaces(surfaces, ray);
    %Now, loop through for each subsequent ray, and stop only when there
    %are no more intersections.
    prevrays{r}=rays{r};
    while ~isempty(intersectionOut)
        %The ray interacts with a surface. Compute this and find the ray
        %that results from the interaction. 
        [newray, power, surfaceNums(r)]=...
            interactAllSurfaces(ray,surfaces,power);
        %Append the intersection to the list
        intersections{r}(:,end+1)=intersectionOut;
        %Does the new ray interact with a surface? We need to know this for
        %the next iteration of the while loop. If the new ray doesn't
        %intersect with a surface, then there will be no more iterations.
        [surfaceNum,intersectionOut,jout]=...
            intersecAllSurfaces(surfaces,newray);
        rays{r}=newray; 
        if ~isempty(intersectionOut)
            ray=newray;
            prevrays=rays;
        end
    end
    if showplot
        plotpath(intersections{r},'r');
    end
    
end
if showplot
    %indicate the final direction of all rays
    for r=1:length(rays)
            plotpath([intersections{r}(:,end),rays{r}.t_ray(50)],'b');
    end
end

