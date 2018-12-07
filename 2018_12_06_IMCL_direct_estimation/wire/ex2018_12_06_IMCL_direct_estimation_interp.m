%%%%%%%%%%%%%%%%%%%%
% �ΏہF���C���^�[�Q�b�g�C���E�����F15.0mm
% �ݒ艹���F1580[m/s]����������
% �œ_�����ʒu�Œ�Fy=0
% ���E�ʒu�F���m
% IMCL������0 %�ɌŒ肷��D
% �g�ʒx���v���t�@�C���ɂ�艹������
% �f�q�Ԏ�M�g���ԍ��͎��M���̐��K�����ݑ��ւ�p����
%%%%%%%%%%%%%%%%%%%%

%% �����ݒ�
clear
dst_path = sprintf('H:/result/2018_12_06_IMCL_direct_estimation/wire');
load("H:/data/kwave/config/t_pos_2board.mat");
load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\rfdata.mat")
load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\medium.mat")
load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\kgrid.mat")
load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\sourse_wave.mat")
v_fat = 1450;%[m/s]
v_muscle = 1580;%[m/s]
t_facing_distance = 0.04;%[m]
[num_sample,num_receiver,num_transmitter] = size(rfdata);
num_echo_receiver = num_transmitter;
reference_point = zeros(1,num_echo_receiver);
reference_point_lowerlimit = zeros(1,num_echo_receiver);%�ώ����]���̂��߂�RF�f�[�^�}�X�L���O�Ɏg��
reference_point_upperlimit = zeros(1,num_echo_receiver);%�ώ����]���̂��߂�RF�f�[�^�}�X�L���O�Ɏg��
distance_from_focal_point_all = zeros(1,num_echo_receiver);
num_depth = (t_pos(2,1)-t_pos(2,101))/kgrid.dx/2 - 3;%'3'�Ƃ���̂́C�ŋߐڋ�����0.4 mm�ł��邱�Ƃ��l�����Ă���D
focal_depth = kgrid.dx*151;
focal_point = [0 ;kgrid.x_vec(300)];
element_pitch = abs(t_pos(1,1) - t_pos(1,2));
lateral_range_max = 0;
lateral_range_min = 0;
num_lateral = round((lateral_range_max -lateral_range_min) / element_pitch)+1;%0.0~4.8 mm�܂�
ind_max_signal = 0;
max_signal = 0;
displacement = zeros(1,num_echo_receiver);
displacement_ref = zeros(1,num_echo_receiver);
error_wavefront = 0;
% focal_signal_total = zeros(1,num_depth);

%% �t�H�[�J�V���O
%�쓮�f�q�̌���
v_reference = v_muscle;
target_element = find((-focal_depth/2<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth/2)));
%��M�p�̎Q�Ɠ_�Z�o
for jj = 1:num_echo_receiver
    distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point);
    delay_time_all = floor(((distance_from_focal_point_all - focal_depth)/v_reference)/kgrid.dt);%[sample]
    reference_point(1,jj) = floor(delay_time_all(1,jj)+(2*focal_depth/v_reference)/kgrid.dt+25-1);
end
%���M�r�[���t�H�[�~���O�i���ʁj
focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
%��M�r�[���t�H�[�J�V���O
focal_signal_total =  receive_focus(focused_rfdata,target_element,reference_point);
interp_rfdata = zeros(num_sample*4,num_echo_receiver);
for jj = 1:num_echo_receiver
    interp_rfdata(:,jj) = interp(focused_rfdata(:,jj),4);
end
source_wave = interp(source_wave,4);

%% ��������
%RF�f�[�^�}�X�L���O(�v�Z�̈�)
min_mask = (min(reference_point) - 10)*4;
max_mask = (min(reference_point) + 100 - 1)*4;
mask_rfdata = zeros(size(interp_rfdata));
mask_rfdata(min_mask:max_mask,:) = 1;
focused_rfdata_mask = mask_rfdata .* interp_rfdata;
min_reference = (min(reference_point) - 30)*4;
max_reference = (min(reference_point) + 500 - 1)*4;
reference_rfdata = zeros(size(interp_rfdata));
reference_rfdata(min_reference:max_reference,:) = 1;
focused_rfdata_reference = reference_rfdata .* interp_rfdata;
% �g�ʌ`�󐄒�
delay_profile = zeros(1,length(target_element));
for jj = 1:length(target_element)-1
    [acor,lag] = xcorr(focused_rfdata_mask(:,target_element(jj+1)),focused_rfdata_reference(:,target_element(jj)),'coeff');
    [~,I] = max(abs(acor));
    displacement(1,target_element(jj)+1) = lag(I);
    delay_profile(1,jj+1) = sum(displacement(1,target_element(1):target_element(jj)+1));
end
ind_central_element = find(target_element == 50);
delay_profile = delay_profile + abs(delay_profile(1,ind_central_element));
[~,delay_offset] = findpeaks(abs(hilbert(focused_rfdata_mask(:,50))),'NPeaks',1,'SortStr','descend');
[~,source_wave_offset] = findpeaks(abs(hilbert(source_wave)),'NPeaks',1,'SortStr','descend');
delay_offset = delay_offset - source_wave_offset;
delay_profile = delay_profile + delay_offset/2;
poly_delay_profile_fitted = polyfit(t_pos(1,target_element),(delay_profile*kgrid.dt/4).^2,2);
displacement_ref(1,target_element(2:end)) = diff(reference_point(target_element));
estimated_velocity = 1/sqrt(poly_delay_profile_fitted(1));
f = polyval(poly_delay_profile_fitted.',t_pos(1,target_element).');
T = table(t_pos(1,target_element).',(delay_profile*kgrid.dt/4).^2.',f,(delay_profile*kgrid.dt/4).^2.'-f,'VariableNames',{'X','Y','Fit','FitError'});
%% �ۑ���
delay_profile = zeros(1,length(target_element));
for jj = 1:length(target_element)-1
    [acor,lag] = xcorr(focused_rfdata_mask(:,target_element(jj+1)),focused_rfdata_reference(:,target_element(jj)),'coeff');
    [~,I] = max(abs(acor));
    displacement(1,target_element(jj)+1) = lag(I);
    delay_profile(1,jj+1) = sum(displacement(1,target_element(1):target_element(jj)+1));
end
delay_profile = delay_profile + abs(delay_profile(1,ind_central_element));
[~,delay_offset] = findpeaks(abs(hilbert(focused_rfdata_mask(:,50))),'NPeaks',1,'SortStr','descend');
[~,source_wave_offset] = findpeaks(abs(hilbert(source_wave)),'NPeaks',1,'SortStr','descend');
delay_profile = delay_profile + delay_offset;
figure;
imagesc(t_pos(1,:)*1e3,kgrid.t_array*1e9,abs(hilbert(interp_rfdata)));
hold on
scatter(t_pos(1,(min(target_element):max(target_element)))*1e3,delay_profile*kgrid.dt/4*1e9,'blue','filled');
scatter(t_pos(1,(min(target_element):max(target_element)))*1e3,(reference_point(1,(min(target_element):max(target_element)))+3)*kgrid.dt*1e9,'red');
hold off
xlabel('lateral[mm]');
ylabel('time[ns]');
axis square;
axis tight;
xlim([t_pos(1,(min(target_element))-1)*1e3 t_pos(1,(max(target_element)+1))*1e3])
ylim([min_reference*kgrid.dt*1e9/4 max_mask*kgrid.dt*1e9/4])
colormap(bone);
colorbar;
caxis([0 max(max(abs(hilbert(focused_rfdata_mask))))])
savefig([dst_path,'\rfdata_detail_interp_rawsignal.fig'])
exportfig([dst_path,'\rfdata_detail_interp_rawsignal'],'png',[400,400])

figure;
imagesc(t_pos(1,:)*1e3,kgrid.t_array*1e9,abs(hilbert(interp_rfdata)));
hold on
scatter(t_pos(1,(min(target_element):max(target_element)))*1e3,delay_profile*kgrid.dt/4*1e9,'blue','filled');
scatter(t_pos(1,(min(target_element):max(target_element)))*1e3,(reference_point(1,(min(target_element):max(target_element)))+3)*kgrid.dt*1e9,'red');
hold off
xlabel('lateral[mm]');
ylabel('time[ns]');
axis square;
axis tight;
xlim([t_pos(1,(min(target_element))-1)*1e3 t_pos(1,(max(target_element)+1))*1e3])
ylim([0 max_mask*kgrid.dt*1e9/2])
colormap(bone);
colorbar;
caxis([0 max(max(abs(hilbert(focused_rfdata_mask))))])
savefig([dst_path,'\rfdata_whole_interp_rawsignal.fig'])
exportfig([dst_path,'\rfdata_whole_interp_rawsignal'],'png',[400,400])
