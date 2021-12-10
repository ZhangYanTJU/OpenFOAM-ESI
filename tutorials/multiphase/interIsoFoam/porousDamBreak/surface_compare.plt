reset

set terminal pngcairo font "helvetica,20" size 50.0in,70.0in

set xrange [0:1]
set yrange [0:0.25]
set xlabel "x [m]"
set ylabel "y [m]"

set grid

graphName="surface_compare.png"
set output graphName

set style rect fc lt -1 fs solid 0.15 noborder
set object  1 rect from 0.29, graph 0 to 0.59, graph 1

array exp_data[10]= ['4', '6', '8', '10', '12', '14', '16', '18', '20', '22']
array time[10]= ['0.4', '0.6', '0.8', '1', '1.2', '1.4', '1.6', '1.8', '2', '2.2']

set multiplot layout 5,2 margins 0.2, 0.8 , 0.2, 0.8 spacing 0.03,0.05

do for[k=1:10] {
set title "t = ".time[k]." [s]"
plot \
    "../porousDamBreak/postProcessing/freeSurface/".time[k]."/freeSurface.raw" ps 3 w points  t "OpenFOAM", \
    "../porousDamBreak/experimentalData/data_".exp_data[k].".dat" w points ps 3 pt 7 t "Liu et. al" 
}
unset multiplot
