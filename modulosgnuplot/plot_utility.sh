#!/bin/bash

help_panel(){
    echo -e "PLOT UTILITY"
    echo -e "\tplot_utility: Herramienta que permite realizar diferentes figuras utilizando gnuplot"
    echo -e "SINOPSIS"
    echo -e "\tplot_utility [-c \"archivo de configuración\"] [-t \"tipo de figura\"]"
    echo -e "OPCIONES"
    echo -e "\t-c \"archivo de configuración\""
    echo -e "\t\tEl archivo de configuración debe ser .txt y tener el siguiente formato:\n"
    echo -e "\t\tPara el tipo [1] se provee un solo archivo de configuración de la forma:\n"
    echo -e "\t\t\tterminal=algún terminal"
    echo -e "\t\t\tsize=[tamaño x, tamaño y]"
    echo -e "\t\t\tdatafile=../direccion/del/archivo.dat"
    echo -e "\t\t\tfont=Fuente, tamaño"
    echo -e "\t\t\toutputfile=Archivo de salida.pdf"
    echo -e "\t\t\tcols=[columnas a graficar]"
    echo -e "\t\t\tlabels=[nombres de labels]"
    echo -e "\t\t\tlegends=[leyendas de cada línea]"
    echo -e "\t\t\tnplots=Número de líneas\n"
    echo -e "\t\tPara el tipo [2] se provee un archivo de texto con los nombres de
                los archivos de configuración:\n"
    echo -e "\t\t\tconfig_1.txt"
    echo -e "\t\t\tconfig_2.txt"
    echo -e "\t\t\tconfig_n.txt\n"
    echo -e "\t\tPara el tipo [3] se debe pasar un archivo similar al del tipo [1], con la diferencia
                de que se tendrán múltiples archivos de datos en un arreglo.\n"
    echo -e "\t-t \"tipo de figura\"\n"
    echo -e "\t\tfigura simple: una figura con una o multiples líneas del mismo archivo de datos [1]"
    echo -e "\t\tfigura multiple: multiples figuras con una o multiples líneas del mismo o múltiples archivos de datos [2]"
    echo -e "\t\tfigura multiple: una figura con una o múltiples líneas de distintos archivos de datos [3]"
}

exec_plot_1(){
    cfile=$1
    type=$2
    arguments=$(grep '=' $cfile | awk -F'=' '{print $2}')
    # gnuplot -c home/juan/cursos/modulosgnuplot/gnup_sm.gp $arguments
    gnuplot -c /home/juan/cursos/modulosgnuplot/gnup_test.gp $arguments
}

exec_plot_2(){
    cfiles=$1
    type=$2
    while read line;
    do
        arguments=$(grep '=' $line | awk -F'=' '{print $2}')
        # gnuplot -c /home/juan/cursos/modulosgnuplot/gnup_sm.gp $arguments
        gnuplot -c /home/juan/cursos/modulosgnuplot/gnup_test.gp $arguments
    done < $cfiles
}

exec_plot_3(){
    cfile=$1
    type=$2
    arguments=$(grep '=' $cfile | awk -F'=' '{print $2}')
    gnuplot -c /home/juan/cursos/modulosgnuplot/gnup_mm.gp $arguments
}


if [ -z $1 ]; then help_panel; exit 1; fi
declare -i parameter_counter=0
while getopts ":c:t:h" arg; do
    case $arg in
        t)  type=$OPTARG
            parameter_counter+=1
            ;;
        c)
            cfile=$OPTARG
            parameter_counter+=1
            ;;
        h)
            help_panel
            exit 0
            ;;
    esac
done

if [ $parameter_counter == 2 ]; then
    if  [ $type == 1 ]; then
        exec_plot_1 $cfile $type
    elif [ $type == 2 ]; then
        exec_plot_2 $cfile $type
    elif [ $type == 3 ]; then
        exec_plot_3 $cfile $type
    fi
else
    help_panel
    exit 1
fi
exit 0