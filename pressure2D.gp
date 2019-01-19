set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
    size 16,9

set output "debug2D0.pdf"

plot[-0.1:1.1][-0.1:1.1] "debug2D/debug2D0" u 1:2:7 with image

set output "debug2D1.pdf"

plot[-0.1:1.1][-0.1:1.1] "debug2D/debug2D1" u 1:2:7 with image

set output "debug2D2.pdf"

plot[-0.1:1.1][-0.1:1.1] "debug2D/debug2D2" u 1:2:7 with image

set output "debug2D845.pdf"

plot[-0.1:1.1][-0.1:1.1] "debug2D/debug2D845" u 1:2:7 with image
