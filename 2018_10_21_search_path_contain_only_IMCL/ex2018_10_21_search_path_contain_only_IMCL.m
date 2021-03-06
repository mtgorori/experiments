%%%%%%%%%%%%%%%%%%%%%%%%%%%%％％％％％％％%%%%%
% ch 50直下のピクセルにフォーカスをかける．
% フォーカス深度を１グリッドずつ変化させて，開口合成受信データを生成する．
% F値を固定する．最近接距離をいくつに設定する？
% x-axis: 20 mmのときに全素子を使うという前提を設ける．
% 対象にする媒質データ：case26を用いる．
% IMCL割合を0 %に固定する．
% kgrid.x_vec が0となるのは251番目の要素．
% グリッド幅が0.1 mm
% 素子間ピッチが0.4 mm
% よって，最近接深度(最近接距離は0.4 mm)
% 深度ごとの細かい素子割当は，floor()を用いる．
% 整相加算の際に参照する音速を1580 m/sとする．
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
load("H:\data\kwave\config\t_pos_2board.mat");
load("H:\data\kwave\result\2018_10_21_highFrequency_variousIMCL\case26_IMCL4.0_5MHz\rfdata.mat");
load("H:\data\kwave\result\2018_10_21_highFrequency_variousIMCL\case26_IMCL4.0_5MHz\medium.mat");
load("H:\data\kwave\result\2018_10_21_highFrequency_variousIMCL\case26_IMCL4.0_5MHz\kgrid.mat");
load("H:\data\kwave\result\2018_10_21_highFrequency_variousIMCL\case26_IMCL4.0_5MHz\sourse_wave.mat");
% load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\rfdata.mat");
% load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\medium.mat");
% load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\kgrid.mat");
% load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\sourse_wave.mat");

% 初期設定
v_fat = 1450;%[m/s]
v_muscle = 1580;%[m/s]
t_facing_distance = 0.04;%[m]
rate_IMCL = [0, 2, 4, 6, 8];
num_depth = (t_pos(2,1)-t_pos(2,101))/kgrid.dx - 3;%'3'とあるのは，最近接距離が0.4 mmであることを考慮している．
[num_sample,num_receiver,num_transmitter] = size(rfdata);
num_echo_receiver = num_transmitter;
reference_point = zeros(2,num_echo_receiver,num_depth);
reference_point_base = linspace(1,100,100)';
reference_point_base = repmat(reference_point_base,1,num_depth);
reference_point(1,:,:) = reference_point_base;
distance_from_focal_point_all = zeros(1,num_echo_receiver);
focal_signal = zeros(length(rate_IMCL),num_depth);
focal_amp = zeros(length(rate_IMCL),num_depth,num_echo_receiver);
focal_phase = zeros(length(rate_IMCL),num_depth,num_echo_receiver);
for kk = 1:length(rate_IMCL)
    v_reference = v_muscle*(1-rate_IMCL(kk)/100) + v_fat*(rate_IMCL(kk)/100);
    focused_rfdata = zeros(num_sample,num_echo_receiver,num_depth);
    focused_rfdata_amp = zeros(num_sample,num_echo_receiver,num_depth);
    focused_rfdata_phase = zeros(num_sample,num_echo_receiver,num_depth);
    for ii = 1:num_depth
        focal_depth = (ii+3)*kgrid.dx;
        focal_point = [0;t_pos(2,1)-focal_depth];
        target_element = find((-focal_depth/2<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth/2)));
        distance_from_focal_point = zeros(1,length(target_element));
        for jj = 1:num_echo_receiver
            distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point);
            delay_time_all = round(((distance_from_focal_point_all - focal_depth)/v_reference)/kgrid.dt);%[sample]
            reference_point(2,jj,ii) = round(delay_time_all(1,jj)+1+(2*focal_depth/v_reference)/kgrid.dt+25);
            %25はfocal_amplitudeを最大にするオフセット．
        end
        for jj = 1:length(target_element)
            distance_from_focal_point(1,jj) = distance_from_focal_point_all(1,target_element(jj));
            % 遅延処理
            delay_time = delay_time_all(1,target_element(jj));%[sample]
            read_range_rfdata = length(delay_time+1:num_sample);
            focused_rfdata(1:read_range_rfdata,:,ii) = focused_rfdata(1:read_range_rfdata,:,ii)...
                +  rfdata(delay_time+1:num_sample,1:100,target_element(jj));%整相加算
        end
        hilb_rfdata = hilbert(focused_rfdata(:,:,ii));
        focused_rfdata_amp(:,:,ii) = abs(hilb_rfdata);
        focused_rfdata_phase(:,:,ii) = atan(imag(hilb_rfdata)./real(hilb_rfdata));
        for jj = 1:length(target_element)
            tmp = focused_rfdata(reference_point(2,target_element(1,jj),ii),target_element(1,jj),ii)/length(target_element);
            focal_signal(kk,ii) = focal_signal(kk,ii)+ tmp;
            focal_amp(kk,ii,target_element(jj)) = focused_rfdata_amp(reference_point(2,target_element(1,jj),ii),target_element(1,jj),ii)/length(target_element);
            focal_phase(kk,ii,target_element(jj)) = focused_rfdata_phase(reference_point(2,target_element(1,jj),ii),target_element(1,jj),ii);
        end
        focal_signal(kk,ii) = abs(hilbert(focal_signal(kk,ii)));
    end
end