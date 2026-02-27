#!/bin/bash


#a=("dataCLASH/datphys_M1931.dat")
a=("dataCLASH/datphys_A209.dat"  \
   "dataCLASH/datphys_A383.dat" \
   "dataCLASH/datphys_M329.dat" \
   "dataCLASH/datphys_M1115.dat" \
   "dataCLASH/datphys_M1931.dat" \
   "dataCLASH/datphys_MS2137.dat" \
   "dataCLASH/datphys_R2129.dat" \
   "dataCLASH/datphys_R2248.dat" \
   "dataCLASH/datphys_M1206.dat")

#b=("parsgT/parsgT_M1931.txt") 
b=("parsgT/parsgT_A209.txt" \
 "parsgT/parsgT_A383.txt" \
   "parsgT/parsgT_M329.txt" \
   "parsgT/parsgT_M1115.txt" \
   "parsgT/parsgT_M1931.txt" \
   "parsgT/parsgT_MS2137.txt" \
   "parsgT/parsgT_R2129.txt" \
   "parsgT/parsgT_R2248.txt" \
   "parsgT/parsgT_M1206.txt")

#c=("OutClash/PM1931Bur.dat")
c=("Output/A209_NFW.dat" \
   "Output/A383_NFW.dat" \
   "Output/M329_NFW.dat" \
   "Output/M1115_NFW.dat" \
   "Output/M1931_NFW.dat" \
   "Output/MS2137_NFW.dat" \
   "Output/R2129_NFW.dat" \
   "Output/R2248_NFW.dat" \
   "Output/M1206_NFW.dat")
   
d=("new_kde/A209_NFW.dat" \
   "new_kde/A383_NFW.dat" \
   "new_kde/M329_NFW.dat" \
   "new_kde/M1115_NFW.dat" \
   "new_kde/M1931_NFW.dat" \
   "new_kde/MS2137_NFW.dat" \
   "new_kde/R2129_NFW.dat" \
   "new_kde/R2248_NFW.dat" \
   "new_kde/M1206_NFW.dat")

# Controllo che gli array abbiano la stessa lunghezza
if [ "${#a[@]}" -ne "${#b[@]}" ] || [ "${#a[@]}" -ne "${#c[@]}" ]; then
    echo "Errore: gli array 'a', 'b' e 'c' devono avere la stessa lunghezza."
    exit 1
fi

for i in "${!a[@]}"; do
    echo "Iterazione $i:"
    echo "  Prima riga: ${a[i]}"
    echo "  Seconda riga: ${b[i]}"
    echo "  Settima riga: ${c[i]}"
	echo "  Undicesima riga: ${d[i]}"
	
	sed -e "1s|.*|${a[i]}|" \
    -e "2s|.*|${b[i]}|" \
    -e "7s|.*|${c[i]}|"\
    -e "11s|.*|${d[i]}|" gomamposst_x.inp > gomamposst_x.tmp && mv gomamposst_x.tmp gomamposst_x.inp

	
    # Qui puoi mettere eventuali comandi da eseguire sul file modificato
    ./sw/bin/gomamposstopt < gomamposst_x.inp
done
