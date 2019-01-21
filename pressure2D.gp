#set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
#    size 10,9

#set output "Velocitymag.pdf"

set term gif animate delay 1

set output "pressure2.gif"

#set cbrange[0:3]

do for [j=0:100] {
  plot[0.35:0.6][0.35:0.6] sprintf("debug2D/debug2D%i",5*j) u 1:2:5 with image title sprintf("%i",5*j)
;}
