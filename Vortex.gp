set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
   size 10,9

set palette defined (0 0.9 0.9 1, 35 0.3 0.3 1, 50 0.6 0.15 0.4, 70 'red', 100 'yellow')

set cbrange[0:1]

set output "ImagesDiapo/VelVortexTheta0PC.pdf"
set title "Vitesse pour {/Symbol q} = 0, 100² mailles, PC"
plot[0:1][0:1] "Vortextheta0PC/it_6747" u 1:2:5 with image

set output "ImagesDiapo/VelVortexTheta1PC.pdf"
set title "Vitesse pour {/Symbol q} = 1, 100² mailles, PC"
plot[0:1][0:1] "Vortextheta1PC/it_6747" u 1:2:5 with image

set output "ImagesDiapo/VelVortexTheta2PC.pdf"
set title "Vitesse pour {/Symbol q} = 2, 100² mailles, PC"
plot[0:1][0:1] "Vortextheta2PC/it_6747" u 1:2:5 with image

set output "ImagesDiapo/VelVortexRusa100.pdf"
set title "Vitesse pour 100² mailles, Rusanov"
plot[0:1][0:1] "VortexRusanov100/it_4356" u 1:2:5 with image

set output "ImagesDiapo/VelVortexRusa200.pdf"
set title "Vitesse pour 200² mailles, Rusanov"
plot[0:1][0:1] "VortexRusanov200/it_8790" u 1:2:5 with image

set output "ImagesDiapo/VelVortexPC200Theta2.pdf"
set title "Vitesse pour {/Symbol q} = 2, 200² mailles, PC"
plot[0:1][0:1] "VortexPC200Theta2/it_13484" u 1:2:5 with image

set cbrange[0:0.026]

set output "ImagesDiapo/MachVortexTheta0PC.pdf"
set title "Nombre de Mach pour {/Symbol q} = 0, 100² mailles, PC"
plot[0:1][0:1] "Vortextheta0PC/it_6747" u 1:2:6 with image

set output "ImagesDiapo/MachVortexTheta1PC.pdf"
set title "Nombre de Mach pour {/Symbol q} = 1, 100² mailles, PC"
plot[0:1][0:1] "Vortextheta1PC/it_6747" u 1:2:6 with image

set output "ImagesDiapo/MachVortexTheta2PC.pdf"
set title "Nombre de Mach pour {/Symbol q} = 2, 100² mailles, PC"
plot[0:1][0:1] "Vortextheta2PC/it_6747" u 1:2:6 with image

set output "ImagesDiapo/MachVortexRusa100.pdf"
set title "Nombre de Mach pour 100² mailles, Rusanov"
plot[0:1][0:1] "VortexRusanov100/it_4356" u 1:2:6 with image

set output "ImagesDiapo/MachVortexRusa200.pdf"
set title "Nombre de Mach pour 200² mailles, Rusanov"
plot[0:1][0:1] "VortexRusanov200/it_8790" u 1:2:6 with image

set output "ImagesDiapo/MachVortexPC200Theta2.pdf"
set title "Nombre de Mach pour {/Symbol q} = 2, 200² mailles, PC"
plot[0:1][0:1] "VortexPC200Theta2/it_13484" u 1:2:6 with image
