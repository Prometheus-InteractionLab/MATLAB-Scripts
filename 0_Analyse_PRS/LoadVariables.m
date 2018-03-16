clear all

prompt = {'number of data sets to load','name of data'};
        dlg_title = 'Number of conditions';
        num_lines = 1;
        answer = inputdlg(prompt,dlg_title,num_lines,{'';''});
        loop = str2num(answer{1});
        name = answer{2};
        
for i = 1:loop
    [data,ivalid] = function_load_processed_data();
end

save(name, '-struct', 'data')

%PlotOTron_V1_6