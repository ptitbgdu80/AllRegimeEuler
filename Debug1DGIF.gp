n=1700
reset
set term gif animate delay 1

set style line 1 lt rgb "#DB0000" lw 2 # m Rouge
set style line 2 lt rgb "#00009E" lw 2 # f Bleu foncé
set style line 3 lt rgb "#00994D" lw 2 # x1 Vert
set style line 4 lt rgb "#F25900" lw 2 # h dashed Orange
set style line 5 lt rgb "#DB006E" lw 2 # j Rouge rose
set style line 6 lt rgb "#7ACAFF" lw 2 # d Bleu clair
set style line 7 lt rgb "#FF66FF" lw 2 # Rose

set output "debug1D.gif"
#set cbrange[25000:250000]

#do for [j=57:n] {
#  plot[0.5:0.54][0.2:0.22] sprintf("debug1D/it_%i",j) u 1:2 w l ls 1 title sprintf("%i",j)
#;}


#do for [j=0:n] {
#  plot[0.3:0.7][0.09:1.01] sprintf("debug1D/it_%i",j) u 1:2 w l ls 1 title sprintf("%i",j)
#;}

do for [j=0:n] {
  plot[0:1][0:310] sprintf("test1D/it_%i",20*j) u 1:4 w l ls 1 title sprintf("%i",20*j)
;}
