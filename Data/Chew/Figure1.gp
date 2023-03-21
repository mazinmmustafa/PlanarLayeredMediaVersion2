reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure1.tex'

set key top left

set xrange [-3:+3]
# set yrange [0:1]

set xtics 1
set ytics 50

# set mxtics 1
# set mytics 1

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

#set arrow from 0.0, 0.0 to 20.0, 0.0 nohead linestyle 1 lc #'black' lw 2

set xlabel '$x$ [m]'
set ylabel '$|E_{x}|$ [V/m]'
set title ''

plot 'DataEJ.dat' using 1:2 with lines lw 2 dt 1 lc 'black' notitle

