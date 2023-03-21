reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure6.tex'

set key bottom left

set xrange [-1000:+1000]
# set yrange [0:1]

# set xtics 10
# set ytics 0.2

# set mxtics 1
# set mytics 1

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

#set arrow from 0.0, 0.0 to 20.0, 0.0 nohead linestyle 1 lc #'black' lw 2

set xlabel '$z$ [nm]'
set ylabel '$\eta_{0}I^{h}_{v}$'
set title ''

plot 'Data.dat' using 1:12 with lines lw 2 dt 1 title '$\mathcal{R}e$', \
     'Data.dat' using 1:13 with lines lw 2 dt 1 title '$\mathcal{I}m$'

