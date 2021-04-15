# gnuplot script to generate polynomial fit

# Set terminal
set terminal png font arial 12 size 800,600

# Dynamic Viscosity
set output 'dynamic-viscosity-fit.png'
# Set plot decorations
set title 'Dynamic Viscosity'
set ylabel 'Dynamic Viscosity [kg/(m-s)]'
set xlabel 'Temperature [K]'
set style line 1 ps 1.5 pt 7 lc 'red'
set style line 2 lw 1.5 lc 'blue'
set grid
set key right top
# Fit data
f(x) = A*x**3 + B*x**2 + C*x + D
fit f(x) 'water-properties.dat' using 2:3 via A, B, C, D
# Save polynomial fit coefficients and plot data
parameter_values = sprintf(" f(x)= %e x^3 + %e x^2\n         %e x + %e", A, B, C, D)
set object 1 rect from 320,0.00158 to 369, 0.00142 fc rgb 'white' 
set label 1 at 320,0.00152 parameter_values
plot 'water-properties.dat' using 2:3 ls 1 t 'Data', f(x) ls 2 t 'Polynomial Fit'

# Thermal Conductivity
set output 'thermal-conductivity-fit.png'
# Set plot decorations
set title 'Thermal Conductivity'
set ylabel 'Thermal Conductivity [W/(m-K)]'
set xlabel 'Temperature [K]'
set style line 1 ps 1.5 pt 7 lc 'red'
set style line 2 lw 1.5 lc 'blue'
set grid
set key left top
# Fit data
f(x) = A*x**3 + B*x**2 + C*x + D
fit f(x) 'water-properties.dat' using 2:4 via A, B, C, D
# Save polynomial fit coefficients and plot data
parameter_values = sprintf(" f(x)= %e x^3 %e x^2\n        + %e x %e", A, B, C, D)
set object 1 rect from 271,0.665 to 318, 0.652 fc rgb 'white' 
set label 1 at 271,0.66 parameter_values
plot 'water-properties.dat' using 2:4 ls 1 t 'Data', f(x) ls 2 t 'Polynomial Fit'

# Specific Heat Capacity
set output 'specific-heat-capacity-fit.png'
# Set plot decorations
set title 'Specific Heat Capacity'
set ylabel 'Specific Heat Capacity [W/(m-K)]'
set xlabel 'Temperature [K]'
set style line 1 ps 1.5 pt 7 lc 'red'
set style line 2 lw 1.5 lc 'blue'
set grid
set key right top
# Fit data
f(x) = A*x**3 + B*x**2 + C*x + D
fit f(x) 'water-properties.dat' using 2:5 via A, B, C, D
# Save polynomial fit coefficients and plot data
parameter_values = sprintf(" f(x)= %e x^3 + %e x^2\n         %e x + %e", A, B, C, D)
set object 1 rect from 315,4214.5 to 362, 4210.5 fc rgb 'white' 
set label 1 at 315,4213 parameter_values
plot 'water-properties.dat' using 2:5 ls 1 t 'Data', f(x) ls 2 t 'Polynomial Fit'

# Density
set output 'density-fit.png'
# Set plot decorations
set title 'Density'
set ylabel 'Density [kg/m^3]'
set xlabel 'Temperature [K]'
set style line 1 ps 1.5 pt 7 lc 'red'
set style line 2 lw 1.5 lc 'blue'
set grid
set key right top
# Fit data
f(x) = A*x**3 + B*x**2 + C*x + D
fit f(x) 'water-properties.dat' using 2:6 via A, B, C, D
# Save polynomial fit coefficients and plot data
parameter_values = sprintf(" f(x)= %e x^3 %e x^2\n         + %e x + %e", A, B, C, D)
set object 1 rect from 322,999.5 to 369,995.5 fc rgb 'white' 
set label 1 at 322,998 parameter_values
plot 'water-properties.dat' using 2:6 ls 1 t 'Data', f(x) ls 2 t 'Polynomial Fit'
