function rays=createrays(numRays,startPoint,centerAngle,fanAngle,mode,...
    totalfanPower,indexOfRefraction)
%Creates a cell array of rays
%numRays: number of rays in the fan
%startPoint: the starting point of the fan
%centerAngle: angle of the middle of the fan [°CCW from x axis]
%fanAngle: opening angle of the fan
%mode: 'linear' creates rays with linearly distributed angles
%      'random' creates ray with random angles in the range


switch mode
    case 'linear'
        angles=linspace(centerAngle-fanAngle/2,...
            centerAngle+fanAngle/2,numRays);
    case 'random'
        angleStart=centerAngle-fanAngle/2;
        angleEnd=centerAngle+fanAngle/2;
        angles=angleStart+rand(1,numRays)*(angleEnd-angleStart); 
    otherwise
        beep
        rays = {};
        return
end

%create an array of direction vectors
raypower=totalfanPower/numRays;
rays=cell(1,numRays);
for i=1:numRays
    direction=[cos(angles(i)*pi/180);...
        sin(angles(i)*pi/180)];
    rays{i}=createray(startPoint,'dir',direction,raypower,...
        indexOfRefraction);
end
