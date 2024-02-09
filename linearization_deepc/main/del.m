data_del = zeros(4,3);
Neq = 1e7;
m  = 10;
plotChannels = {'RotSpeed','PtfmPitch','GenPwr','TwrBsMyt','RotTorq','RootMyc1'};

figure

for LC = 1:4

    switch LC
        case 1
            outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs3_tp12\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
            file1 = load('fullSimFinal2\steadyWind_irrWaves_Hs3Tp12_Mp_preview_QP_scaled.mat');
        case 2
            outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs4p3_tp10\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
            file1 = load('fullSimFinal2\steadyWind_irrWaves_Hs4p3Tp10_Mp_preview_QP_scaled.mat');
        case 3
            outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs3_tp12_turbwind\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
            file1 = load('fullSimFinal2\turbWind_irrWaves_Hs3Tp12_Fsg_Mp_preview_QP_scaled.mat');
        case 4
            outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs4p3_tp10_turbwind\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
            file1 = load('fullSimFinal2\turbWind_irrWaves_Hs4p3Tp10_Fsg_Mp_preview_QP_scaled.mat');
    end

    
    [data, channels, units, headers] = ReadFASTtext(outFile);
    id = find(ismember(channels,plotChannels{4}));
    twr = data(:,id);

    id = find(ismember(channels,plotChannels{5}));
    lss = data(:,id);

    id = find(ismember(channels,plotChannels{6}));
    oop = data(:,id);

    twr_deepc = file1.out(4,:).*file1.Dy(4,4);
    lss_deepc = file1.out(5,:).*file1.Dy(5,5);
    oop_deepc = file1.out(6,:).*file1.Dy(6,6);

    plot(twr_deepc)
    hold on


    S = max(twr)-min(twr);
    del_twr_baseline = (S^m/Neq)^(1/m);

    S = max(twr_deepc)-min(twr_deepc)
    del_twr_deepc = (S^m/Neq)^(1/m)

    data_del(LC,1) = del_twr_deepc/del_twr_baseline - 1;

    S = max(lss)-min(lss);
    del_lss_baseline = (S^m/Neq)^(1/m);

    S = max(lss_deepc)-min(lss_deepc);
    del_lss_deepc = (S^m/Neq)^(1/m);

    data_del(LC,2) = del_lss_deepc/del_lss_baseline - 1;

    S = max(oop)-min(oop);
    del_oop_baseline = (S^m/Neq)^(1/m);

    S = max(oop_deepc)-min(oop_deepc);
    del_oop_deepc = (S^m/Neq)^(1/m);

    data_del(LC,3) = del_oop_deepc/del_oop_baseline - 1;


end

X = categorical({'Tower base loads','LSS torque','OoP blade moment'});
X = reordercats(X,{'Tower base loads','LSS torque','OoP blade moment'});

figure
bar(X,data_del)
legend('Load case 1', 'Load case 2', 'Load case 3', 'Load case 4','Location','northwest')
ylabel('Performance ratio (-)')
grid on
set(gcf,'Color','White')



