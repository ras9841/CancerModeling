#!bin/bash
#Configs[1]="identical.cfg" 
#Tags[1]="identical" 
Configs[1]="diff_prop.cfg" 
Tags[1]="diff_prop" 
Configs[2]="diff_elast.cfg" 
Tags[2]="diff_elast" 
Configs[3]="diff_surf_E.cfg" 
Tags[3]="diff_surf_E" 
#Configs[5]="diff_bulk_mod.cfg" 
#Tags[5]="diff_bulk_mod" 
Configs[4]="full.cfg" 
Tags[4]="full" 
Configs[5]="full_div.cfg"
Tags[5]="division"

NUM=5
TESTS=3

for ((i=1; i<=$NUM; i++ ))
do
    sh run_many.sh ${Configs[i]} ${Tags[i]} ${TESTS}
done
