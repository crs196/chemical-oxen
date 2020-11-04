function discharge_coeff=cd_chart_calc(LD,CD,lift,diameter)
non_dim_lift=round(lift/diameter,4);
discharge_coeff=CD(LD==non_dim_lift);
if isempty(discharge_coeff)
    discharge_coeff=0;
end