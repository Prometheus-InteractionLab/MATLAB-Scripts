function data_p = function_volts_2_phys(data_a,x,y,yo)

% Transducer sensitivity [(phys. unit)/V], form:  sens = polyval(pcylS,Tpcylxd);
sens = polyval(x,y);

% Transducer offset [V], form:  offset = polyval(pcylO,Tpcylxd);
offset = polyval(yo,y);

% Apply transducer calibration, form:  pcylA = sens.*(pcylA-offset);
data_p = sens.*(data_a-offset);
