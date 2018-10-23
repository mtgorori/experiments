%�쐬�F2018/09/26-----------------------------------------------------------------------------------
%�ŒZ�o�H���̊e�f�q�y�A�ɂ�����X�̌o�H�Ɋւ��āC�e�}���̐�L�ʐϔ䂩��Z�o����
%���ω�����z��v(���M�f�q�ԍ�,�}���ԍ�,���g��)�Ɋi�[����X�N���v�g�D
%------------------------------------------------------------------------------------------------------
clear;close all
load("H:/data/kwave/result/2018_09_08_variousFrequency/statistics.mat",'tof_cell')
load("H:/data/kwave/param/param_2board.mat")
%% �����ݒ�
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
%% lt(=v)�f�[�^�̍쐬
v = zeros(t_num,model_num,freq_num);
for ii = 1:freq_num
    for jj = 1:model_num
        for kk = 1:t_num
            v(kk,jj,ii) = leng(100+kk,kk)/tof_cell(100+kk,kk,jj,ii);
        end
    end
end
figure;
plot(v(:,1,2));
figure;
waterfall(v(:,:,1)');
xlabel('�f�q�ԍ�')
ylabel('�}���ԍ�')
zlabel('���ω���[m/s]')
exportfig("H:\result\2018_09_26_ltplot\lt_from_signal_waterfall",'png',[400,400])
figure;
imagesc(v(:,:,1)');
xlabel('�f�q�ԍ�')
ylabel('�}���ԍ�')
axis tight
axis square
c = colorbar;
c.Label.String = '���ω���[m/s]';
exportfig("H:\result\2018_09_26_ltplot\lt_from_signal_imagesc",'png',[300,300])
%% �Q�Ɨp��lt(=v)�f�[�^�̍쐬
v_reference = zeros(t_num,model_num);
% �K�؂ȎQ�Ɣ͈͂�T���D
displace_t_m = abs(t_pos(1,2) - t_pos(1,1));%[m]
displace_t = displace_t_m / param.grid.dx;%[t_num]
% ��f�q�ɂ��O���b�h2���̉�f�l�𕽋ω����ĎQ�Ɨp��lt(=v)�f�[�^�̍쐬���s���D
load("H:/data/kwave/medium/2018_08_10_realisticScatter/case1.mat")
figure;
imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,medium.sound_speed);
ax = gca;
ax.YDir = 'normal';
xlabel('x��[mm]')
ylabel('y��[mm]')
exportfig("H:\result\2018_09_26_ltplot\medium_case1",'png',[300,300]);
%  �f�q���ǂ̃O���b�h�Ɉʒu���Ă��邩�����߂�D
ind_t_pos = zeros(2,200);
[~,ind_t_pos(1,:)] = min(abs(t_pos(1,:)-kgrid.y_vec));
[~,ind_t_pos(2,:)] = min(abs(t_pos(2,:)-kgrid.x_vec));
ref_medium = medium;
%�ǂ̃O���b�h���Q�Ƃ��ĎQ�Ɨp��lt�f�[�^���쐬���Ă��邩���m�F����D�܂��C���s�I�Ɉ�Ԗڂ̃��f���Ɋւ��ĎQ�Ɨp�̃f�[�^���쐬����D
for i = 1:t_num
    v_reference(i,1) = mean2(ref_medium.sound_speed(ind_t_pos(1,i)-1:ind_t_pos(1,i)+1,ind_t_pos(2,100+i):ind_t_pos(2,i)));
    ref_medium.sound_speed(ind_t_pos(1,i)-1:ind_t_pos(1,i)+1,ind_t_pos(2,100+i):ind_t_pos(2,i)) = 0;
end
figure;
imagesc(kgrid.x_vec*1000,kgrid.y_vec*1000,ref_medium.sound_speed);
ax = gca;
ax.YDir = 'normal';
xlabel('x��[mm]')
ylabel('y��[mm]')
exportfig("H:\result\2018_09_26_ltplot\refer_range_medium_case1",'png',[300,300]);
v_reference = zeros(t_num,model_num);
% ��f�l���ω����s��
for jj = 1:model_num
    loadfilename = sprintf("H:/data/kwave/medium/2018_08_10_realisticScatter/case%d.mat",jj);
    load(loadfilename);
    ref_medium = medium;
    for ii = 1:t_num
        v_reference(ii,jj) = mean2(ref_medium.sound_speed(ind_t_pos(1,ii)-1:ind_t_pos(1,ii)+1,ind_t_pos(2,100+ii):ind_t_pos(2,ii)));
    end
end
v_reference = repmat(v_reference,1,1,6);
figure;
waterfall(v_reference(:,:,1)');
xlabel('�f�q�ԍ�')
ylabel('�}���ԍ�')
exportfig("H:\result\2018_09_26_ltplot\lt_from_medium_waterfall_2MHz",'png',[300,300]);
figure;
imagesc(v_reference(:,:,1)');
xlabel('�f�q�ԍ�')
ylabel('�}���ԍ�')
axis tight
axis square
c = colorbar;
c.Label.String = '���ω���[m/s]';
exportfig("H:\result\2018_09_26_ltplot\lt_from_medium_imagesc_2MHz",'png',[300,300])
%% v�v���b�g�̔�r
%water-fall�v���b�g
figure;
subplot(1,2,1);
waterfall(v_reference(:,:,1)');
xlabel('�f�q�ԍ�')
ylabel('�}���ԍ�')
zlabel('���ω���[m/s]')
title({'�f�q�y�A���Ƃ̕��ω����l';'(�}���f�[�^����쐬)'})
axis square
subplot(1,2,2);
waterfall(v(:,:,1)');
xlabel('�f�q�ԍ�')
ylabel('�}���ԍ�')
zlabel('���ω���[m/s]')
title({'�f�q�y�A���Ƃ̕��ω����l';'(RF�f�[�^����쐬)'})
axis square
exportfig("H:\result\2018_09_26_ltplot\lt_compared_waterfall_2MHz",'png',[800,400])
%imagesc�v���b�g
figure;
subplot(2,2,1);
imagesc(v_reference(:,:,1)');
xlabel('�f�q�ԍ�')
ylabel('�}���ԍ�')
title({'�f�q�y�A���Ƃ̕��ω����l';'(�}���f�[�^����쐬)'})
axis square
colorbar
caxis([1450 1580])
ax = gca;
ax.YDir = 'normal';
subplot(2,2,2);
imagesc(v(:,:,1)');
xlabel('�f�q�ԍ�')
ylabel('�}���ԍ�')
title({'�f�q�y�A���Ƃ̕��ω����l';'(RF�f�[�^����쐬)'})
axis square
colorbar
caxis([1450 1580])
ax = gca;
ax.YDir = 'normal';
subplot(2,2,[3,4]);
imagesc(abs(v(:,:,1)'-v_reference(:,:,1)'));
xlabel('�f�q�ԍ�')
ylabel('�}���ԍ�')
title('�f�q�y�A���Ƃ̕��ω����l����덷')
axis square
colorbar
ax = gca;
ax.YDir = 'normal';
exportfig("H:\result\2018_09_26_ltplot\lt_compared_imagesc_2MHz",'png',[800,400])

frq = [2000, 1000, 500, 200, 100, 50];
figure;
for i = 1:6
    subplot(2,3,i)
    imagesc(v(:,:,i)'-v_reference(:,:,i)');
    xlabel('�f�q�ԍ�')
    ylabel('�}���ԍ�')
    title('�o�H���ω�������덷')
    axis square
    colorbar
    ax = gca;
    ax.YDir = 'normal';
    caxis([-35 35])
    mytitlename = sprintf('���S���g��=%d kHz',frq(i));
    title({'�o�H���ω�������덷';mytitlename});
end
exportfig("H:\result\2018_09_26_ltplot\lt_compared_imagesc_multi_freq",'png',[900,500])
%% �e��f�[�^�ۑ�
