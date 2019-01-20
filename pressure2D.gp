#set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
#    size 10,9

#set output "Velocitymag.pdf"

set term gif animate delay 100

set output "pressure2.gif"

set cbrange[0:0.4]

do for [j=0:2] {
  plot[-0.1:1.1][-0.1:1.1] sprintf("debug2D/debug2D%i",j) u 1:2:4 with image title sprintf("%i",j)
;}
