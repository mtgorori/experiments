clear;close all
load("//Azlab-fs01/��������/�lwork/�|��(��)/data/kwave/result/2018_09_08_variousFrequency/statistics.mat",'tof_cell')
load("//Azlab-fs01/��������/�lwork/�|��(��)/data/kwave/param/param_2board.mat")
load("//Azlab-fs01/��������/�lwork/�|��(��)/data/kwave/result/2018_09_08_variousFrequency/statistics.mat",'frq')

t_size = 40.e-3;
t_num = 200;%�g�����X�f���[�T��
t_pos = zeros(2, t_num);%�Z���T�ʒu
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%�f�q������������[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;
[r_num,t_num,model_num,freq_num] = size(tof_cell);
leng = zeros(r_num,t_num);
for i = 1:r_num
    for j = 1:t_num
        leng(i,j) = sqrt((t_pos(1,i)-t_pos(1,j))^2+(t_pos(2,i)-t_pos(2,j))^2);
    end
end


v = zeros(t_num,model_num,freq_num);
for ii = 1:freq_num
    for jj = 1:model_num
        for kk = 1:t_num
            v(kk,jj,ii) = leng(100+kk,kk)/tof_cell(100+kk,kk,jj,ii);
        end
    end
end

load("//Azlab-fs01/��������/�lwork/�|��(��)/data/kwave/medium/2018_08_10_realisticScatter/case1.mat")
v_reference = zeros(t_num,model_num);
ind_t_pos = zeros(2,200);
[~,ind_t_pos(1,:)] = min(abs(t_pos(1,:)-kgrid.y_vec));
[~,ind_t_pos(2,:)] = min(abs(t_pos(2,:)-kgrid.x_vec));
ref_medium = medium;
% ��f�l���ω����s��
for jj = 1:model_num
    loadfilename = sprintf("//Azlab-fs01/��������/�lwork/�|��(��)/data/kwave/medium/2018_08_10_realisticScatter/case%d.mat",jj);
    load(loadfilename);
    ref_medium = medium;
    for ii = 1:t_num
        v_reference(ii,jj) = mean(ref_medium.sound_speed(ind_t_pos(1,ii),ind_t_pos(2,100+ii):ind_t_pos(2,ii)));
    end
end
v_reference = repmat(v_reference,1,1,6);
%imagesc�v���b�g
titlename = ["2 MHz","1 MHz", "500 kHz", "200 kHz", "100 kHz", "50 kHz"];
for ii = 1:6
    figure;
    subplot(1,3,1);
    imagesc((((v_reference(:,:,ii)-1580)/(1450-1580))*100)');
    xlabel('�f�q�ԍ�')
    ylabel('�}���ԍ�')
    title('�f�q�y�A���Ƃ̐���IMAT�ܗL��')
    axis square
    colorbar
    c = colorbar;
    c.Label.String = '[%]';
    caxis([0 100])
    ax = gca;
    ax.YDir = 'normal';
    subplot(1,3,2);
    imagesc((((v(:,:,ii)-1580)/(1450-1580))*100)');
    xlabel('�f�q�ԍ�')
    ylabel('�}���ԍ�')
    title('�f�q�y�A���Ƃ̐���IMAT�ܗL��')
    axis square
    colorbar
    c = colorbar;
    c.Label.String = '[%]';
    caxis([0 100])
    ax = gca;
    ax.YDir = 'normal';
    subplot(1,3,3);
    imagesc(abs((((v(:,:,ii)-1580)/(1450-1580))*100)'-(((v_reference(:,:,ii)-1580)/(1450-1580))*100)'));
    xlabel('�f�q�ԍ�')
    ylabel('�}���ԍ�')
    title({'�f�q�y�A���Ƃ�';'IMAT�ܗL������덷'})
    axis square
    colorbar
    c = colorbar;
    c.Label.String = '[%]';
    caxis([0 50])
    ax = gca;
    ax.YDir = 'normal';
    axes;
    title(titlename(ii),'FontSize',18,'FontName', 'Cambria');
    axis off;
    
    figname= sprintf("H:/result/2018_09_26_ltplot/lt_compared_imagesc_%d",ii);
    
    exportfig(figname,'png',[1000,300])
end