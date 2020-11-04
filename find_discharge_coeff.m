%%RUN THIS FILE

valve=input('Intake (1) or exhaust (2)');
switch valve
    case 1
        [LD,CD]=cd_chart_intake;
    case 2
        [LD,CD]=cd_chart_exhaust;
end

for i=1:length(lift)
discharge(i,1)=cd_chart_calc(LD,CD,lift(i),1.1875);
end