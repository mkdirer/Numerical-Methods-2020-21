set terminal post enhanced colour solid font 25  # wybor formatu, w jakim chcemy utworzyc wynikowy rysunek
set encoding iso_8859_2

set style line 1 pt 7 lc rgb "blue"          
set style line 2 pt 7 lc rgb "goldenrod"  
set style line 3 pt 7 lc rgb "dark-spring-green"               
set style line 4 pt 7 lc rgb "red" 
set style line 5 pt 7 lc rgb "turquoise"              
set style line 6 pt 7 lc rgb "brown" 

# results.txt -- plik z pierwiastkami wartosci wlasnych. Zawiera 7 kolumn: parametr alfa oraz 6 odpowiadajacych mu pierwiastkow wartosci wlasnych.
set output "chart2.eps" # Pierwiastki wartosci wlasnych w funkcji parametru alpha
set xlabel "{/Symbol a}"
set ylabel "{/Symbol w}({/Symbol a})"
set key outside horizontal top center
set logscale x
set xtics ("2" 2, "10" 10, "100" 100)
#set multiplot
plot for [k=0:5] "results.txt" u 1:(column(k+2)) w lp ls k+1 t "{/Symbol w}_".k.""
# Wstawiony rysunek wewnetrzny:
