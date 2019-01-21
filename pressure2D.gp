#set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
#    size 10,9

#set output "Velocitymag.pdf"

set term gif animate delay 5

set output "test1D.gif"

set cbrange[10000:100000]

do for [j=0:1000] {
  plot[0.3:0.7][-0.01:1.01] sprintf("test1D/it_%i",20*j) u 1:2:4 with image title sprintf("%i",20*j)
;}
