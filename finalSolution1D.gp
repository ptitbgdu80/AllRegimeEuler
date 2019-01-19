set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
    size 16,9

set output "image.pdf"

plot[-0.01:1.01][9050:101000] "resultats/resultats520" u 1:6 w l
