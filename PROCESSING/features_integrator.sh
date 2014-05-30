#! bin/sh

# CAMBIAR FILTROS Y LA CARPETA DE SALIDA!!

for file in GENE.SHORTINTRON.65_8-30 #EXON.SHORTINTRON.65_8-30
	do
		awk -F'[\t]' '{ if (a[$2]++ == 0) print }' ../processed_data/features/${file}/${file} > ${file}_shortintrons #tiene sentido a nivel de exon cuando cogemos los short introns como neutros. Ojo que si cogemos los 4 fold a nivel de exon esto implica que solo cogemos la info del ultimo exon!!! Es para evitar que multiplique los intrones cortos tantas veces como exones tiene un gen. A nivel de gen no importa coger el short introns o el 4 fold. 
		for filters in `less G31` #change de FILTERSS!!!!!! --> E2 (CHR:MUTATION:RECOMBINATION:TRANSCRIPTS)
			do
				echo ${filters}
				perl features_integrator_${file}.pl -coding ../processed_data/features/${file}/${file} -shortintrons ${file}_shortintrons -filter ${filters} > ./LIST44/${filters} #change the neutral class site!!!				
			done
		rm ${file}_shortintrons
	done

exit 0

## awk -F'[\t]' '{ if (a[$2]++ == 0) print }' ../../INTEGRATION/features/${file} > ${file}_shortintrons