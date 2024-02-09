% for idxSim = 3:4
% 
% switch idxSim
%     case 4
%         baselineController1
%     case 3
%         baselineController2
%     case 2
%         main_turbWind1
%     case 1
%         main_turbWind2
% end
% 
% end

for idxSim = 1:3

switch idxSim
    case 1
        getOLdata
    case 2
        main_steadyWind1
    case 3
        main_steadyWind2
end

end