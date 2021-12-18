clear all


result_folder='D:\CST Studio\Comparisons_Hook_Coupler_Designs\Eigenmode_C3794_Qext_Optimized_Coupler_4hook_from_450MHz\Export';



c0=299792458;
mu0=4*pi*1e-7;
imp_cal=1;
freq_all=importdata(strcat(result_folder,'\Frequency (Multiple Modes).txt'));

Mode_num=size(freq_all,1);
freq_all=freq_all(:,end)*1e6;%second column is frequency: 1e6 for MHz and 1e9 for GHz

for i=1:Mode_num
Ex{1,i}=importdata(strcat(result_folder,sprintf('\\e_X (Z)_Mode %u.txt',i)));
Ey{1,i}=importdata(strcat(result_folder,sprintf('\\e_Y (Z)_Mode %u.txt',i)));
Hx{1,i}=importdata(strcat(result_folder,sprintf('\\h_X (Z)_Mode %u.txt',i)));
Hy{1,i}=importdata(strcat(result_folder,sprintf('\\h_Y (Z)_Mode %u.txt',i)));
Z_axis=Ex{1,i}(:,1)*1e-3;%1e-3 to convert mm to m
VT_x(i,1)=trapz(Z_axis,(Ex{1,i}(:,2)+1i*Ex{1,i}(:,3)-c0*mu0*(Hy{1,i}(:,2)+1i*Hy{1,i}(:,3))).*exp(1i*2*pi*freq_all(i)*Z_axis/c0));
VT_y(i,1)=trapz(Z_axis,(Ey{1,i}(:,2)+1i*Ey{1,i}(:,3)+c0*mu0*(Hx{1,i}(:,2)+1i*Hx{1,i}(:,3))).*exp(1i*2*pi*freq_all(i)*Z_axis/c0));

end

R_Q_x=abs(VT_x.^2)./(2*pi*freq_all);
R_Q_y=abs(VT_y.^2)./(2*pi*freq_all);
R_Q_T=R_Q_x+R_Q_y;

if imp_cal==1
  Qext_all=importdata(strcat(result_folder,'\Qext.txt'));
  Qext_all=Qext_all(:,2);
  Z_T=Qext_all.*R_Q_T*0.5*2*pi.*freq_all/c0;
  Z_export=[freq_all R_Q_T Qext_all Z_T];
  
  Z_x=Qext_all.*R_Q_x*0.5*2*pi.*freq_all/c0;
  Z_y=Qext_all.*R_Q_y*0.5*2*pi.*freq_all/c0;

  
end

fs = linspace(470e6, 600e6, 2000);%the frequency interval to calculate impedance from eigenmodes

for i=1:Mode_num
    Z_sum_x(i,:)=1*Z_x(i)./(1+1i*Qext_all(i)*(fs/freq_all(i)-freq_all(i)./fs));
    Z_sum_y(i,:)=1*Z_y(i)./(1+1i*Qext_all(i)*(fs/freq_all(i)-freq_all(i)./fs));

end
Z_x_all=sum(Z_sum_x);%impedance in x direction
Z_y_all=sum(Z_sum_y);%impedance y direction

figure
semilogy(fs,abs(Z_x_all),'-b')
hold on
semilogy(fs,abs(Z_y_all),'-r')
semilogy(freq_all,abs(Z_x),'ob')
semilogy(freq_all,abs(Z_y),'or')


%writematrix(Z_export,'A_trans_secondDipole.csv');
%writematrix(Z_export,'A_trans_FirstDipole.csv');