#!/usr/bin/gnuplot

!rm 'iter'*'.dat'
!cat log | grep "fInf =" | awk '{print sqrt($7*$7)}' > 'iter_fInf.dat'
!cat log | grep "fInf =" | awk '{print $13}' > 'iter_Uinf.dat'
!cat log | grep "fInf =" | awk '{print sqrt($45*$45)}' > 'iter_FtotArea.dat'
!cat log | grep "fInf =" | awk '{print sqrt($51*$51)}' > 'iter_FtotVol.dat'


Re=0.0156905
Mg=0.0145282
Pr=0.462963
uinf=0.104603
uinf_corr = 1*uinf

# Re=7.84524
# Mg=7.26411
# Pr=0.462963
# uinf=52.3016
# uinf_corr = 0.82855189982776*uinf

# Re=78.4524
# Mg=72.6411
# Pr=0.462963
# uinf=523.016
# uinf_corr = 0.575*uinf

set title sprintf("2D: Stationary bubble terminal velocity: Re = %g, Mg = %g, Pr = %g, U_{inf} = %g mm/s", Re, Mg, Pr, uinf)

set xlabel "Accumulated PIMPLE iterations"
set ylabel "Force / (10^{-12} N)"
set y2label "Velocity / (mm/s)"

set xrange [0:*]
# set xrange [0:80000]
# set yrange [*:100]
# set y2range [0:0.14]

set ytics nomirror
set ytics
set y2tics

set format y "10^{%L}"
set logscale y


set terminal pdf size 18 cm, 14 cm
set output 'plot.pdf'

plot 0 w l ls -1 notitle, \
     uinf w l lw 2 lc black dt 2 title "analytical" axes x1y2, \
     uinf_corr w l lw 2 lc black dt 4 title "analytical, corrected" axes x1y2, \
     'iter_FtotArea.dat' u ($1*1e+9) w l lw 2 t "F_{tot, A} / 1000" axes x1y1, \
     'iter_FtotVol.dat' u ($1*1e+9) w l lw 2 dt 3 t "F_{tot, V} / 1000" axes x1y1, \
     'iter_fInf.dat' u ($1*1e+12) w l lw 2 t "F_{inf}" axes x1y1, \
     'iter_Uinf.dat' u ($1*-1e+3) w l lw 2 t "U_{inf}" axes x1y2