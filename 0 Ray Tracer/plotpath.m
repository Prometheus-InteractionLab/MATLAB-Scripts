function total_length=plotpath(pointlist,color)
%makes a new figure and displays the points given in pointlist
%color: rgb vector or text (i.e. 'r' or [1 0 0])
%pointlist: list of points to plot in which each point is a two-element 
%           column vector [x;y]
%total_length: the total distance between the points

plot(pointlist(1,:),pointlist(2,:),'Color',color,'LineWidth',1,...
    'Marker','o','MarkerSize',2)
axis equal
set(gca,'FontSize',12);
set(gca,'FontWeight','bold');
xlabel('x','FontWeight','bold','FontSize',16);
ylabel('y','FontWeight','bold','FontSize',16);



%calculate total path length from start to end. define function reference
%to determine length between two points
% pathlength=@ (x1,x2) sqrt((x2(1,1)-x1(1,1))^2+(x2(2,1)-x1(2,1))^2);
% %loop and calculate distances between pairs of points
% total_length=0;
% for ii=1:size(pointlist,2)-1
%     total_length=total_length+...
%         pathlength(pointlist(:,ii),pointlist(:,ii+1));
% end
