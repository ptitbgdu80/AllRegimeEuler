set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
   size 10,9

set palette defined (0 0.9 0.9 1, 35 0.3 0.3 1, 50 0.6 0.15 0.4, 70 'red', 100 'yellow')

set style line 1 lt rgb "#DB0000" lw 2 # m Rouge
set style line 2 lt rgb "#00009E" lw 2 # f Bleu foncé
set style line 3 lt rgb "#00994D" lw 2 # x1 Vert
set style line 7 lt rgb "#F25900" lw 2 # h dashed Orange
set style line 5 lt rgb "#DB006E" lw 2 # j Rouge rose
set style line 6 lt rgb "#7ACAFF" lw 2 # d Bleu clair
set style line 4 lt rgb "#FF66FF" lw 2 # Rose

set cbrange[0:1.75]

set output "ImagesDiapo/VelRiemann2Dtheta0.pdf"
set title "Vitesse pour {/Symbol q} = 0, 200² mailles, PC"
plot[0:1][0:1] "Riemann2Dtheta0/it_3283" u 1:2:5 with image

set output "ImagesDiapo/VelRiemann2Dtheta1.pdf"
set title "Vitesse pour {/Symbol q} = 1, 200² mailles, PC"
plot[0:1][0:1] "Riemann2Dtheta1/it_3252" u 1:2:5 with image

set output "ImagesDiapo/VelRiemann2Dtheta2.pdf"
set title "Vitesse pour {/Symbol q} = 2, 200² mailles, PC"
plot[0:1][0:1] "Riemann2Dtheta2/it_3251" u 1:2:5 with image

set output "ImagesDiapo/VelRiemann2DRusa.pdf"
set title "Vitesse pour 200² mailles, Rusanov"
plot[0:1][0:1] "Riemann2DRusanov/it_853" u 1:2:5 with image

set cbrange[0:3.15]

set output "ImagesDiapo/MachRiemann2Dtheta0.pdf"
set title "Nombre de Mach pour {/Symbol q} = 0, 200² mailles, PC"
plot[0:1][0:1] "Riemann2Dtheta0/it_3283" u 1:2:6 with image

set output "ImagesDiapo/MachRiemann2Dtheta1.pdf"
set title "Nombre de Mach pour {/Symbol q} = 1, 200² mailles, PC"
plot[0:1][0:1] "Riemann2Dtheta1/it_3252" u 1:2:6 with image

set output "ImagesDiapo/MachRiemann2Dtheta2.pdf"
set title "Nombre de Mach pour {/Symbol q} = 2, 200² mailles, PC"
plot[0:1][0:1] "Riemann2Dtheta2/it_3251" u 1:2:6 with image

set output "ImagesDiapo/MachRiemann2DRusa.pdf"
set title "Nombre de Mach pour 200² mailles, Rusanov"
plot[0:1][0:1] "Riemann2DRusanov/it_853" u 1:2:6 with image



set output "ImagesDiapo/Density_y=x.pdf"
set title "Densité en y=x"
plot[-0.01:1.01][0.1:1.6] "Riemann2Dtheta0/y=x" u 1:2 w l ls 1 t "{/Symbol q} = 0", "Riemann2Dtheta1/y=x" u 1:2 w l ls 2 t "{/Symbol q} = 1", "Riemann2Dtheta2/y=x" u 1:2 w l ls 3 t "{/Symbol q} = O(M)", "Riemann2DRusanov/y=x" u 1:2 w l ls 4 t "Rusanov"

set output "ImagesDiapo/Pressure_y=x.pdf"
set title "Pression en y=x"
plot[-0.01:1.01][0:1.6] "Riemann2Dtheta0/y=x" u 1:3 w l ls 1 t "{/Symbol q} = 0", "Riemann2Dtheta1/y=x" u 1:3 w l ls 2 t "{/Symbol q} = 1", "Riemann2Dtheta2/y=x" u 1:3 w l ls 3 t "{/Symbol q} = O(M)", "Riemann2DRusanov/y=x" u 1:3 w l ls 4 t "Rusanov"

set output "ImagesDiapo/Velocity_Mag_y=x.pdf"
set title "Vitesse en y=x"
plot[-0.01:1.01][0:1.8] "Riemann2Dtheta0/y=x" u 1:4 w l ls 1 t "{/Symbol q} = 0", "Riemann2Dtheta1/y=x" u 1:4 w l ls 2 t "{/Symbol q} = 1", "Riemann2Dtheta2/y=x" u 1:4 w l ls 3 t "{/Symbol q} = O(M)", "Riemann2DRusanov/y=x" u 1:4 w l ls 4 t "Rusanov"

set output "ImagesDiapo/Mach_Number_y=x.pdf"
set title "Nombre de Mach en y=x"
plot[-0.01:1.01][0:3.2] "Riemann2Dtheta0/y=x" u 1:5 w l ls 1 t "{/Symbol q} = 0", "Riemann2Dtheta1/y=x" u 1:5 w l ls 2 t "{/Symbol q} = 1", "Riemann2Dtheta2/y=x" u 1:5 w l ls 3 t "{/Symbol q} = O(M)", "Riemann2DRusanov/y=x" u 1:5 w l ls 4 t "Rusanov"


set output "ImagesDiapo/Density_x=0.75.pdf"
set title "Densité en x=0.75"
plot[-0.01:1.01][0.5:1.55] "Riemann2Dtheta0/x=0.75" u 1:2 w l ls 1 t "{/Symbol q} = 0", "Riemann2Dtheta1/x=0.75" u 1:2 w l ls 2 t "{/Symbol q} = 1", "Riemann2Dtheta2/x=0.75" u 1:2 w l ls 3 t "{/Symbol q} = O(M)", "Riemann2DRusanov/x=0.75" u 1:2 w l ls 4 t "Rusanov"

set output "ImagesDiapo/Pressure_x=0.75.pdf"
set title "Pression en x=0.75"
plot[-0.01:1.01][0.2:1.6] "Riemann2Dtheta0/x=0.75" u 1:3 w l ls 1 t "{/Symbol q} = 0", "Riemann2Dtheta1/x=0.75" u 1:3 w l ls 2 t "{/Symbol q} = 1", "Riemann2Dtheta2/x=0.75" u 1:3 w l ls 3 t "{/Symbol q} = O(M)", "Riemann2DRusanov/x=0.75" u 1:3 w l ls 4 t "Rusanov"

set output "ImagesDiapo/Velocity_Mag_x=0.75.pdf"
set title "Vitesse en x=0.75"
plot[-0.01:1.01][0:1.3] "Riemann2Dtheta0/x=0.75" u 1:4 w l ls 1 t "{/Symbol q} = 0", "Riemann2Dtheta1/x=0.75" u 1:4 w l ls 2 t "{/Symbol q} = 1", "Riemann2Dtheta2/x=0.75" u 1:4 w l ls 3 t "{/Symbol q} = O(M)", "Riemann2DRusanov/x=0.75" u 1:4 w l ls 4 t "Rusanov"

set output "ImagesDiapo/Mach_Number_x=0.75.pdf"
set title "Nombre de Mach en x=0.75"
plot[-0.01:1.01][0:1.4] "Riemann2Dtheta0/x=0.75" u 1:5 w l ls 1 t "{/Symbol q} = 0", "Riemann2Dtheta1/x=0.75" u 1:5 w l ls 2 t "{/Symbol q} = 1", "Riemann2Dtheta2/x=0.75" u 1:5 w l ls 3 t "{/Symbol q} = O(M)", "Riemann2DRusanov/x=0.75" u 1:5 w l ls 4 t "Rusanov"
