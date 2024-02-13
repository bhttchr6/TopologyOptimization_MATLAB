function T = stimulus_function(STIMULUS_PARAM, t)
Tmax = STIMULUS_PARAM.Tmax;
Tmin = STIMULUS_PARAM.Tmin;
tinitial = STIMULUS_PARAM.tinitial;
HR = STIMULUS_PARAM.HR;
CR = STIMULUS_PARAM.CR;
tfinal = STIMULUS_PARAM.t_final;
tmax_1 = (Tmax-Tmin)/HR + tinitial;
tmax_2 = tfinal - (Tmax-Tmin)/CR;
Theating = (Tmin + HR*(t-tinitial))*HEAVY(tmax_1-t);
Tcooling = (Tmin - CR*(t-tfinal))*HEAVY(t-tmax_2);

T = Theating + Tmax*2*HEAVY(t-tmax_1) + Tcooling;


end

function val = HEAVY(parameter)
if parameter>0
    val = 1;
end
if parameter <0
    val = 0;
end
if parameter ==0
    val = 0.5;
end
end