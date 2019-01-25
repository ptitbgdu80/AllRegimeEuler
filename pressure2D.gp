set palette defined (0 0.9 0.9 1, 35 0.3 0.3 1, 50 0.6 0.15 0.4, 70 'red', 100 'yellow')

set term gif animate delay 1

set output "VitesseRiemann2DRusanov.gif"

set cbrange[0:1.75]

do for [j=0:1000] {
  plot[0:1][0:1.1] sprintf("Riemann2DRusanov/it_%i",50*j) u 1:2:5 with image title sprintf("%i",50*j)
;}
