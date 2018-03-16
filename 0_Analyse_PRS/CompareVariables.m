fig = figure(4)

 prompt = {'number of data sets to compare','Matlab marker symbol [sdv^xo+*<>ph]','Matlab colors [bgrcmykw], default is random color','plot legend data','round to nearest','plot legend text','legend units'};
        dlg_title = 'Number of conditions';
        num_lines = 1;
        answer = inputdlg(prompt,dlg_title,num_lines,{'';'';'';'Ttdc';'10';'TDC Temp';'K'});
        loop = str2num(answer{1});
        mark = answer{2};
        color = answer{3};
        plotlegend = answer{4};
        digits = str2num(answer{5});
        legendtext = answer{6};
        units = answer{7};
        
for i = 1:loop
    fig,set(gcf,'color','w','name','Comparing Data','Position',[85 160 1456 754]);
    [fig, data] = function_plot_processed_data_AHRRvsgIMEP(fig,mark(i),color(i),plotlegend,legendtext,digits,units);
end