result_folder = 'D:\Users\Shahnam\TwoCellCavity_forSWELL\Qext_WG\SWELL_QEXT_IgorwithTransition\Export';

result_folder = 'D:\Users\Shahnam\TwoCellCavity_forSWELL\Qext_WG\SWELL_QEXT_IgorwithTransition_OpenBC2\Export';
result_folder = 'D:\Users\Shahnam\TwoCellCavity_forSWELL\Qext_WG\SWELL_QEXT_IgorwithTransition_OpenBC\Export';

# result_folder = 'D:\Users\Shahnam\TwoCellCavity_forSWELL\Qext_WG\SWELL_QEXT_NoDielectric_improved\Export';
# result_folder = 'D:\Users\Shahnam\TwoCellCavity_forSWELL\Qext_WG\SWELL_QEXT_NoDielectric_improved2\Export';

# result_folder = 'D:\Users\Shahnam\TwoCellCavity_forSWELL\Qext_WG\SWELL_QEXT_short\Export';
# result_folder = 'D:\Users\Shahnam\TwoCellCavity_forSWELL\Qext_WG\SWELL_QEXT_short_OpenBC2\Export';

c0 = 299792458;
mu0 = 4 * pi * 1e-7;
imp_cal = 1;
freq_all = importdata(strcat(result_folder, '\Frequency (Multiple Modes).txt'));

Mode_num = size(freq_all, 1);
freq_all = freq_all(:, end)*1e6; # second
column is frequency: 1e6
for MHz and 1e9 for GHz

for i=1:Mode_num
Ex
{1, i} = importdata(strcat(result_folder, sprintf('\\e_X (Z)_Mode #u.txt', i)));
Ey
{1, i} = importdata(strcat(result_folder, sprintf('\\e_Y (Z)_Mode #u.txt', i)));
Hx
{1, i} = importdata(strcat(result_folder, sprintf('\\h_X (Z)_Mode #u.txt', i)));
Hy
{1, i} = importdata(strcat(result_folder, sprintf('\\h_Y (Z)_Mode #u.txt', i)));
Z_axis = Ex
{1, i}(:, 1)*1e-3; # 1e-3
to
convert
mm
to
m
VT_x(i, 1) = trapz(Z_axis, (Ex{1, i}(:, 2) + 1
i * Ex
{1, i}(:, 3)-c0 * mu0 * (Hy{1, i}(:, 2) + 1
i * Hy
{1, i}(:, 3))).*exp(1
i * 2 * pi * freq_all(i) * Z_axis / c0));
VT_y(i, 1) = trapz(Z_axis, (Ey{1, i}(:, 2) + 1
i * Ey
{1, i}(:, 3)+c0 * mu0 * (Hx{1, i}(:, 2) + 1
i * Hx
{1, i}(:, 3))).*exp(1
i * 2 * pi * freq_all(i) * Z_axis / c0));

end

R_Q_x = abs(VT_x. ^ 2). / (2 * pi * freq_all);
R_Q_y = abs(VT_y. ^ 2). / (2 * pi * freq_all);
R_Q_T = R_Q_x + R_Q_y;

if imp_cal == 1
    Qext_all = importdata(strcat(result_folder, '\Q-Factor (lossy E) (Multiple Modes).txt'));
    Qext_all = Qext_all(:, 2);
    Z_T = Qext_all. * R_Q_T * 0.5 * 2 * pi. * freq_all / c0;
    Z_export = [freq_all R_Q_T Qext_all Z_T];

end

# writematrix(Z_export, 'C_trans_secondDipole.csv');
# writematrix(Z_export, 'C_trans_FirstDipole.csv');