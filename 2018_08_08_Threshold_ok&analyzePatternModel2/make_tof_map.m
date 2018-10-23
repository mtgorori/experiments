clear; close all;
addpath(genpath('H:\codes'))
cd('H:\data\kwave\param')
load('param_2board.mat')
tof_cell  = zeros(200,100,20,10);%{Receiver, Transmitter, ERate, Num, }
aveSOS = zeros(20,10);%{ERate, Num}
steSOS = zeros(20,10);
leng = zeros(1,100);
t_size = 40.e-3;
t_num = 200;%トランスデューサ数
t_pos = zeros(2, t_num);%センサ位置
t_pos(1,1:t_num/2) = -t_size/2:t_size/(t_num/2-1):t_size/2 ;%素子水平方向距離[m]
t_pos(2,1:t_num/2) = t_size/2;
t_pos(1,t_num/2+1:t_num) = t_pos(1,1:t_num/2);
t_pos(2,t_num/2+1:t_num) = -t_size/2;

for i = 1:10
    for j = 2:2:20
        tic
        cd('H:\data\kwave\result\2018_08_05_patternScatter2')
        myfilename = sprintf('rEMCL%d & nEMCL%d',j,i);
        cd(myfilename)
        load('rfdata.mat')
        load('kgrid.mat')
        tof_data = threshold_picker(rfdata,kgrid);
        cd('H:\result\2018_08_08_analyzePatternModel2');
        save(myfilename,'tof_data');
        tof_cell(:,:,j,i) = tof_data;
        [Min,ind]= min(tof_cell(101:end,:,j,i));
        ind = ind+100;
        for k = 1: t_num/2
            leng(1,k) = norm(t_pos(:,k)-t_pos(:,ind(1,k)));
        end
        aveSOS(j,i) =sum(leng./Min)/(t_num/2);
        steSOS(j,i) =std(leng./Min,0,2)/sqrt(t_num/2);
        toc
    end
end

cd('H:\result\2018_08_08_analyzePatternModel2');
myfilename = sprintf('2018_08_08_aveSOS&steSOS');
save(myfilename,'aveSOS','steSOS','tof_cell');