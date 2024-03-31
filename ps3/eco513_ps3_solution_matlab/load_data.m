% Load from data files
dat_monthly = table2timetable(readtable('eco513_ps3_Monthly.txt'));
dat_quarterly = table2timetable(readtable('eco513_ps3_Quarterly.txt'));
dat = innerjoin(retime(dat_monthly,'quarterly','mean'), dat_quarterly);

% Read data series
price = dat.GDPDEF;
unrate = dat.UNRATE;
infl = [NaN; 100*log(price(2:end)./price(1:end-1))];

% Enforce common sample for data series
data_matrix = [infl unrate];
sample = find(all(~isnan(data_matrix),2),1):find(all(~isnan(data_matrix),2),1,'last');
data_matrix_sample = data_matrix(sample,:);