#!/bin/bash

help_panel(){
    echo -e "PLOT UTILITY"
    echo -e "\tplot_utility: Herramienta que permite realizar diferentes figuras utilizando gnuplot"
    echo -e "SINOPSIS"
    echo -e "\tplot_utility [-t archivo de configuración]"
    echo -e "OPCIONES"
    echo -e "\t-t archivo de configuración"
    echo -e "\t\tEl archivo de configuración debe ser .txt y tener el siguiente formato:"
}

exec_plot(){
    tfile=$1
    arguments=$(grep '=' $tfile | awk -F'=' '{print $2}')
    # echo $arguments
    gnuplot -c ../../../../modulosgnuplot/gnup_sm.gp $arguments
}


if [ -z $1 ]; then help_panel; exit 1; fi
declare -i parameter_counter=0
while getopts ":t:h" arg; do
    case $arg in
        t)
            tfile=$OPTARG
            parameter_counter+=1
            ;;
        h)
            help_panel
            exit 0
            ;;
    esac
done

if [ $parameter_counter == 1 ]; then
    exec_plot $tfile
else
    help_panel
    exit 1
fi
exit 0