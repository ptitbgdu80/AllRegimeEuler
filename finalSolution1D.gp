set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
    size 14,9

set style line 1 lt rgb "#DB0000" lw 2 # m Rouge
set style line 2 lt rgb "#00009E" lw 2 # f Bleu foncé
set style line 3 lt rgb "#00994D" lw 2 # x1 Vert
set style line 4 lt rgb "#F25900" lw 2 # h dashed Orange
set style line 5 lt rgb "#DB006E" lw 2 # j Rouge rose
set style line 6 lt rgb "#7ACAFF" lw 2 # d Bleu clair
set style line 7 lt rgb "#FF66FF" lw 2 # Rose

set output "Density1D.pdf"

plot[-0.01:1.01][0.1:1.01] "Riemann1D_theta0/it_506" u 1:2 w l ls 1 t "{/Symbol q} = 0", "Riemann1D_theta1/it_436" u 1:2 w l ls 2 t "{/Symbol q} = 1", "Riemann1D_theta2/it_469" u 1:2 w l ls 3 t "{/Symbol q} = O(M)"

set output "Pressure1D.pdf"

plot[-0.01:1.01][9050:101000] "Riemann1D_theta0/it_506" u 1:3 w l ls 1 t "{/Symbol q} = 0", "Riemann1D_theta1/it_436" u 1:3 w l ls 2 t "{/Symbol q} = 1", "Riemann1D_theta2/it_469" u 1:3 w l ls 3 t "{/Symbol q} = O(M)"

set output "Velocity_Mag1D.pdf"

plot[-0.01:1.01][-0.1:300] "Riemann1D_theta0/it_506" u 1:4 w l ls 1 t "{/Symbol q} = 0", "Riemann1D_theta1/it_436" u 1:4 w l ls 2 t "{/Symbol q} = 1", "Riemann1D_theta2/it_469" u 1:4 w l ls 3 t "{/Symbol q} = O(M)"

set output "Mach_Number1D.pdf"

plot[-0.01:1.01][-0.01:0.95] "Riemann1D_theta0/it_506" u 1:5 w l ls 1 t "{/Symbol q} = 0", "Riemann1D_theta1/it_436" u 1:5 w l ls 2 t "{/Symbol q} = 1", "Riemann1D_theta2/it_469" u 1:5 w l ls 3 t "{/Symbol q} = O(M)"
