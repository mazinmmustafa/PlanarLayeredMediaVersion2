reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure.tex'

set key bottom right

set xrange [0:90]
set yrange [0:1]

set xtics 10
set ytics 0.2

set mxtics 1
set mytics 1

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

#set arrow from 0.0, 0.0 to 20.0, 0.0 nohead linestyle 1 lc #'black' lw 2

set xlabel '$\theta_{i}$ [deg]'
set ylabel '$|\Gamma|^2$'
set title ''

plot 'Gamma.dat' using 1:2 with lines lw 2 dt 1 lc 'black' title '$e$', \
     'Gamma.dat' using 1:3 with lines lw 2 dt 2 lc 'black' title '$h$'

