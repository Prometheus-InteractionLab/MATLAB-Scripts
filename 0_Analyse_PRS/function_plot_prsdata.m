function function_plot_prsdata(ca,data,cyc,inewfig,varargin)

    set(0,'defaultAxesFontName','Arial')
    set(0,'defaultTextFontName','Arial')

if inewfig
    figure,hold on
end

if nargin<5
cm = colormap('jet');
else
cm = varargin{1};
end

for ii=1:length(cyc)
    icolor = round(63*(ii-1)/(length(cyc)-1)+1);
    plot(ca(cyc(ii),:),data(cyc(ii),:),'o','color',cm(icolor,:),'linestyle','-','linewidth',0.5,'Markersize',1.5,'MarkerFaceColor',cm(icolor,:));
end
   
set(gcf,'color','w')
set(gca,'box','on','linewidth',1,'fontsize',12)
