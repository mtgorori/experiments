% %% 1: 測定した最大音速と理論的な最大音速，そして各素子での最大音速をもつ経路の変動がわかるようなアニメーションを作る．
% %% さらに，最大音速を持つ経路を伝搬する波の強度を表示させて，波の強度が最大音速に対する補正値となりうるかを検討するようなアニメーションにする．
% clear;close all
% load("H:/data/kwave/result/2018_10_15_realisticScatter_variousIMCL_Correct/case1_IMCL0.2/kgrid.mat")
% load("H:/data/kwave/config/t_pos_2board.mat")
% load("H:/data/kwave/result/2018_10_15_realisticScatter_variousIMCL_Correct/EMCL_num_75/statistics.mat")
% [num_IMCL, num_EMCL] = size(rate_EMCLs);
% v_muscle = 1580;%筋肉の音速[m/s]
% v_fat = 1450;%脂肪の音速[m/s]
% reference_max_SOS = zeros(1,10);
% sum_amplitude = zeros(1,75);
% for i = 1:10
%     reference_max_SOS(1,i) = v_muscle*((100-2*i/10)/100) + v_fat*((2*i/10)/100);
% end
% for ii = 1:num_EMCL
%     loadfilename = sprintf("H:/data/kwave/result/2018_10_15_realisticScatter_variousIMCL_Correct/case%d_IMCL%0.1f/rfdata.mat",ii,0.2);
%     load(loadfilename);
%     sum_amplitude(1,ii) = sum(abs(rfdata(:,100+ind_fastest(2,1,ii),ind_fastest(1,1,ii))));
%     disp(ii);
% end
clear mov;
fr(1:num_EMCL*num_IMCL) = struct('cdata',[],'colormap',[]);
for ii = 1:num_EMCL
    for jj = 1:num_IMCL
        figure(1);
        subplot(2,2,1);
        loadfilename = sprintf('H:/data/kwave/result/2018_10_15_realisticScatter_variousIMCL_Correct/case%d_IMCL%0.1f/medium.mat',ii,rate_IMCLs(jj,1)/10);
        load(loadfilename);
        imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
        axis equal
        axis tight
        hold on
        scatter(t_pos(2,:)*1000,t_pos(1,:)*1000,'rs','filled');
        xlabel('y-axis[mm]')
        ylabel('x-axis[mm]')
        hold on
%         for i = 1:length(ind_fastest_all)
%             plot([t_pos(2,i)*1000;t_pos(2,100+ind_fastest_all(i,jj,ii))*1000],[t_pos(1,i)*1000;t_pos(1,100+ind_fastest_all(i,jj,ii))*1000],'.-k','LineWidth',0.001);
%         end
        plot([t_pos(2,ind_fastest(1,jj,ii))*1000;t_pos(2,100+ind_fastest(2,jj,ii))*1000],[t_pos(1,ind_fastest(1,jj,ii))*1000;t_pos(1,100+ind_fastest(2,jj,ii))*1000],'.-r','LineWidth',3)
        hold off
        titlename = sprintf('case%d: EMCL=%0.2f %%',ii, rate_EMCLs(1,ii));
        title(titlename);
        subplot(2,2,2);
        plot(rate_EMCLs(1,:));
        hold on
        scatter(ii,rate_EMCLs(1,ii),'filled','r');
        hold off
        xlabel('medium number');
        ylabel('EMCL[%]')
        subplot(2,2,3);
        plot(sum_amplitude);
        hold on
        scatter(ii, sum_amplitude(1,ii),'filled','red');
        hold off
        xlabel('medium number');
        ylabel('total amplitude');
        subplot(2,2,4);
        plot(0.2:0.2:2.0,max_SOS_all(:,ii));
        hold on
        plot(0.2:0.2:2.0,reference_max_SOS(1,:));
        scatter(jj*0.2, max_SOS_all(jj,ii),'filled','red');
        hold off
        xlabel('IMCL[%]');
        ylabel('maximum velocity[m/s]')
        xlim([0 2.2]);
        ylim([1558 1587])
        legend('measure','reference','Location','best');
        pause(0.05);
        drawnow;
        fr(jj+(ii-1)*num_IMCL) = getframe(gcf);
    end
end
% play movie
figure;
movie(fr)
% save movie
cd('H:/result/2018_10_19_maxSOS_consider_attenuation')
mv = VideoWriter('maximum_velocity_path_IMCLchg_attenuation','MPEG-4');
mv.FrameRate = 8; % ← fpsと同じ %ART:3
open(mv)
writeVideo(mv,fr)
close(mv)