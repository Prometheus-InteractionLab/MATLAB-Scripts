function PID_Estimate()
% Ian Liu, 8/8/2015
% Estimates the engine temperature based on the PID set temperature

%Obtains Mass flow rate from user and rounds it to 30, 40, or 50 g/s (where
%data was taken)
format short g
i=1;
while i~=3
    choice=input('What would you like to do? \n<E>stimate temperature\n<V>iew data\n<Q>uit program\n', 's');
    choice=lower(choice);
    switch choice
        case 'e'
        MassFlow=round(input('What is the mass flow rate (will round to 30,40, or 50 g/s)? ')*.1)/.1;
        if MassFlow <30
         MassFlow = 30;
        elseif MassFlow >50
            MassFlow = 50;
        end %if
        %Column specifies with excel column has the correct data based on mass flow
        column=MassFlow/10-1;
        

        wantedTemp=input('What is the desired air temperature? ');
        plotdata=xlsread('E:\Data\PID Data.xlsx','Summary', 'A3:D20');

        %Determines what two set temperatures the desired will fall between
        lowIndexValue=find(plotdata(:,column)<=wantedTemp);
        if isempty(lowIndexValue) == 1
         lowIndexValue = 1;
        end %if
        lowIndexValue=max(lowIndexValue);
        highIndexValue=lowIndexValue+1;

        %Linearly interpolates the data and calculates the set temperature
        slope=(plotdata(highIndexValue,column)-plotdata(lowIndexValue,column))/(plotdata(highIndexValue,1)-plotdata(lowIndexValue,1));
        setTemp=round(((wantedTemp-plotdata(lowIndexValue,column))/slope+plotdata(lowIndexValue,1))*10)/10

        %Plots correct set of data based on mass flow
        if exist('fig') == 1
            if ishandle(fig)==1
                close(fig);
            end
        end
        fig=figure();
        if MassFlow <=30
         plot(plotdata(:,1),plotdata(:,2),'-ok')
        elseif MassFlow == 40
         plot(plotdata(:,1),plotdata(:,3),'-ok')
        elseif MassFlow >= 50
         plot(plotdata(:,1),plotdata(:,4),'-ok')
        end %if

        %Creates final touches on the plot
        hold on
        x=plot(setTemp,wantedTemp,'^b');
        set(x, 'MarkerSize', 8);
        x=text(setTemp-20,wantedTemp+10, strcat('(',num2str(setTemp),', ', num2str(wantedTemp),')'));
        set(x, 'Color', 'blue')
        legend('Data Points', 'Desired Temperature','Location','southeast')
        legend('boxoff')
        xlabel('Set Temperature (celsius)')
        ylabel('Actual Temperature (celsius)')
        assignin('base', 'Data', plotdata);

        case 'v'
        [num,txt,rawdata]=xlsread('E:\Data\PID Data.xlsx','Summary', 'A1:H8');
        x=size(rawdata);
        fig1=figure('name', 'Raw data');
        plotdata=cell2mat(rawdata(3:x,1:4));
        subplot(2,2,1)
        plot(plotdata(:,1),plotdata(:,2),'-ok')
        title('30 g/s')
        xlabel('Set Temperature (celsius)')
        ylabel('Actual Temperature (celsius)')
        subplot(2,2,2)
        plot(plotdata(:,1),plotdata(:,3),'--sb')
        title('40 g/s')
        xlabel('Set Temperature (celsius)')
        ylabel('Actual Temperature (celsius)')
        subplot(2,2,3)
        plot(plotdata(:,1),plotdata(:,4),':^r')
        title('50 g/s')
        xlabel('Set Temperature (celsius)')
        ylabel('Actual Temperature (celsius)')
        subplot(2,2,4)
        title('Overlayed Lines')
        plot(plotdata(:,1),plotdata(:,3),'--sb',plotdata(:,1),plotdata(:,2),'-ok', plotdata(:,1),plotdata(:,4),':^r')
        xlabel('Set Temperature (celsius)')
        ylabel('Actual Temperature (celsius)')
        legend('30 g/s', '40 g/s', '50 g/s', 'Location', 'northwest')
        legend('boxoff')

        i=1;
        while i==1 
        
        [X Y]=ginput(1);
            if X > 40 && X < 60 && Y > 20 && Y < 55
                plot_temp=xlsread('E:\Data\PID Data.xlsx','Data by Temperature', 'A1:F177');
                mainAir=1;
            elseif X > 80 && X < 120 && Y > 60 && Y < 80
                plot_temp=xlsread('E:\Data\PID Data.xlsx','Data by Temperature', 'G1:M298');
                mainAir=2;
            elseif X > 130 && X < 170 && Y > 85 && Y < 105
                plot_temp=xlsread('E:\Data\PID Data.xlsx','Data by Temperature', 'N1:T308');
                mainAir=3;
            elseif X > 180 && X < 220 && Y > 110 && Y < 125
                plot_temp=xlsread('E:\Data\PID Data.xlsx','Data by Temperature', 'U1:Z164');
                mainAir=4;
            elseif X > 230 && X < 270 && Y > 135 && Y < 155
                plot_temp=xlsread('E:\Data\PID Data.xlsx','Data by Temperature', 'AA1:AE211');
                mainAir=5;
            elseif X > 280 && X < 320 && Y > 165 && Y < 180
                plot_temp=xlsread('E:\Data\PID Data.xlsx','Data by Temperature', 'AF1:AI233');
                mainAir=6;
            else 
                disp('Please choose a point on a graph')
            end %if
            fig2=figure('name', strcat(num2str(cell2mat(rawdata(mainAir+2,1))), 'C'));
            subplot(2,2,1);
            [a b]=size(plot_temp);
            ii=2;
            while ii <= b-2
                hold on
                plot(plot_temp(:,1), plot_temp(:,ii), '-k');
                ii=ii+1
            end
            axis([0 length(plot_temp)/2+20 -inf max(max(plot_temp(:,2:4)))+15]);
            title(strcat('30 g/s' ));
            hold on
            x=[0 300];
            y=[cell2mat(rawdata(mainAir+2,2)) cell2mat(rawdata(mainAir+2,2))];
            plot(x,y);
            xlabel('Time Elapsed (min)')
            ylabel('Actual Temperature (celsius)')
            subplot(2,2,2);
            plot(plot_temp(:,1), plot_temp(:,b-1), '--b');
            axis([0 length(plot_temp)/2+20 -inf max(max(plot_temp(:,2:4)))+15]);
            title('40 g/s')
            hold on
            x=[0 300];
            y=[cell2mat(rawdata(mainAir+2,3)) cell2mat(rawdata(mainAir+2,3))];
            plot(x,y);
            xlabel('Time Elapsed (min)')
            ylabel('Actual Temperature (celsius)')
            subplot(2,2,3)
            plot(plot_temp(:,1), plot_temp(:,b), ':r');
            axis([0 length(plot_temp)/2+20 -inf max(max(plot_temp(:,2:4)))+15]);
            title('50 g/s')
            hold on
            x=[0 300];
            y=[cell2mat(rawdata(mainAir+2,4)) cell2mat(rawdata(mainAir+2,4))];
            plot(x,y);
            xlabel('Time Elapsed (min)')
            ylabel('Actual Temperature (celsius)')
            i=4;
            while i==4
                check=lower(input('<C>hoose another point\n<R>eturn to main menu\n<V>iew detailed data\n<Q>uit program\n', 's'));
                if check == 'c'
                    i=1;
                    if ishandle(fig2)
                        close(fig2);
                    end
                elseif check =='r'
                    if ishandle(fig2)
                       close(fig2);
                    end
                    if ishandle(fig1)
                       close(fig1);
                    end
                    i=2;
                elseif check =='q'
                    i=3;          
                elseif check =='v'
                    [num2,txt2,rawdata2]=xlsread('E:\Data\PID Data.xlsx','Summary 2', 'A1:G53');
                    rawdata2
                    i=4;
                end %if
            end%while
        end %while
        case 'q'
            i=3;
    end %switch
end %while
end %function