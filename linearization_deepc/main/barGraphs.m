%% Clean environment
clearvars;clc;close all;
rng('default')

%%
load('fullSimFinal2\postProc.mat')

data = zeros(size(fieldnames(LC1),1),4);

fn = fieldnames(LC1);
for k=1:numel(fn)
    data(k,1) = (LC1.(fn{k})(2)/LC1.(fn{k})(3) - 1)*1; 
end

fn = fieldnames(LC2);
for k=1:numel(fn)
    data(k,2) = (LC2.(fn{k})(2)/LC2.(fn{k})(3) - 1)*1; 
end

fn = fieldnames(LC3);
for k=1:numel(fn)
    data(k,3) = (LC3.(fn{k})(2)/LC3.(fn{k})(3) - 1)*1; 
end

fn = fieldnames(LC4);
for k=1:numel(fn)
    data(k,4) = (LC4.(fn{k})(2)/LC4.(fn{k})(3) - 1)*1; 
end
%%
X = categorical({'Rotor speed variance','Rotor speed RMSE','Blade pitch variance','Blade pitch ADC','Platform pitch variance','Power variance'});
X = reordercats(X,{'Power variance','Rotor speed variance','Rotor speed RMSE','Blade pitch variance','Blade pitch ADC','Platform pitch variance'});

figure
bar(X,data)
legend('Load case 1', 'Load case 2', 'Load case 3', 'Load case 4','Location','northwest')
ylabel('Performance ratio (-)')
grid on
set(gcf,'Color','White')

% %%
% data_abs = zeros(size(fieldnames(LC1),1),4);
% 
% fn = fieldnames(LC1);
% for k=1:numel(fn)
%     data_abs(k,1) = (LC1.(fn{k})(2) - P0_steady.(fn{k}))/(LC1.(fn{k})(3) - P0_steady.(fn{k})); 
% end
% 
% fn = fieldnames(LC2);
% for k=1:numel(fn)
%     data_abs(k,2) = (LC2.(fn{k})(2) - P0_steady.(fn{k}))/(LC2.(fn{k})(3) - P0_steady.(fn{k})); 
% end
% 
% fn = fieldnames(LC3);
% for k=1:numel(fn)
%     data_abs(k,3) = (LC3.(fn{k})(2) - P0_turb.(fn{k}))/(LC3.(fn{k})(3) - P0_turb.(fn{k})); 
% end
% 
% fn = fieldnames(LC4);
% for k=1:numel(fn)
%     data_abs(k,4) = (LC4.(fn{k})(2) - P0_turb.(fn{k}))/(LC4.(fn{k})(3) - P0_turb.(fn{k})); 
% end
% 
% %%
% X = categorical({'Rotor STD','Rotor RMSE','Blade pitch STD','Blade pitch ADC','Platform pitch STD','Mean power'});
% X = reordercats(X,{'Mean power','Rotor STD','Rotor RMSE','Blade pitch STD','Blade pitch ADC','Platform pitch STD'});
% 
% %% no wind
% figure
% bar(X,data_abs(:,1:2))
% legend('Load case 1', 'Load case 2', 'Load case 3', 'Load case 4','Location','northwest')
% ylabel('Normalised (-)')
% yline(1)
% grid on
% set(gcf,'Color','White')
% 
% %% turb wind
% figure
% bar(X,data_abs(:,3:4))
% legend('Load case 1', 'Load case 2', 'Load case 3', 'Load case 4','Location','northwest')
% ylabel('Normalised (-)')
% yline(1)
% grid on
% set(gcf,'Color','White')

