cd "/Users/ryanbroll/Documents/GitHub/class/eco513/ps5/"

import delimited RNUSBIS.csv, clear
gen date2 = date(date,"YMD")
format date2 %td
gen date3 = mofd(date2)
format date3 %tm
tsset date3

dfuller rnusbis, lags(1) regress
