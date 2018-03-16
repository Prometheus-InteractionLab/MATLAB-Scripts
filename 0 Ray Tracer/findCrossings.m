function [X,Y] = findCrossings(rays)
%takes a cell of rays and computes the points at which they intersect each
%other.
disp('Calculating image locations...')
X=zeros(1,length(rays)^2)/0; %Initialize as vector of NaNs
Y=zeros(1,length(rays)^2)/0; %Initialize as vector of NaNs
cnt=0;
%Find at the crossing of each ray with every other ray
for r1=1:length(rays)
    %Output to the command window every 100 rays
    if mod(r1,100)==0 || r1==1
        fprintf('Ray %0.0f / %0.0f\n',r1,length(rays));
    end
    for r2=1:length(rays)
        if ~(r1==r2)
            ray1=rays{r1};
            ray2=rays{r2};
            m1=ray1.pu(2)/ray1.pu(1);
            m2=ray2.pu(2)/ray2.pu(1);

            x01=ray1.start(1);
            x02=ray2.start(1);
            y01=ray1.start(2);
            y02=ray2.start(2);

            x=( m1*x01 - m2*x02 - (y01-y02) ) / (m1-m2);
            if ~isnan(x)
                cnt=cnt+1;
                X(1,cnt)=x;
                Y(1,cnt)=m2*(x-x02)+y02;
            end
        end
    end
end
X=X(~isnan(X)); %if no crossing is found, remove it
Y=Y(~isnan(Y)); 
disp('Done.')