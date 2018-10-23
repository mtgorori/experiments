cd 'H:\data\kwave\result\2018_08_28_realisticScatter\case1'
load('sourse_wave')
load('param')
load('kgrid')

q_value = zeros(400,100);
for ii = 1:400
    cd 'H:\data\kwave\result\2018_08_28_realisticScatter'
    myfilename = sprintf('case%d',ii);
    cd(myfilename)
    load('rfdata.mat');
    figure;
    disp(ii)
    for jj = 1:100
        %figure保存用
        fs = param.sensor.freq;
        t = kgrid.t_array;
        y = fft(rfdata(:,100+jj,jj));
        n = length(t);
        y0 = fftshift(y);         % shift y values
        f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
        power0 = abs(y0).^2/n;    % 0-centered power
%         plot(f0,(power0))
%         xlabel('Frequency')
%         ylabel('Power')
%         titlename = sprintf('case%d: No.%d',ii,jj);
%         title(titlename);
%         cd 'H:\result\2018_09_04_Qvalue_400samples'
%         savefilename1 = sprintf('H:/result/2018_09_04_Qvalue_400samples/case%d_no%d',ii,jj);
%         exportfig(savefilename1,'png',[200,200]);
        %q値算出用
        [~,ind_maxf0] = max(power0((n+1)/2+1:end));
        f0_max = f0(ind_maxf0+(n+1)/2);
        f_width = FWHM(f0((n+1)/2+1:end),power0((n+1)/2+1:end));
        q_value(ii,jj) = f0_max/f_width;
        disp(jj)
    end
end
cd 'H:\result\2018_09_04_Qvalue_400samples'
save('q_value.mat','q_value')
