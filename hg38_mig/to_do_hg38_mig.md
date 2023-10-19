- ask David
	- we are NOT using the mask!!!
		- Jesus: "Respecto a la máscara, creo que no hace falta usarla. La máscara está hecha para evitar los sitios con bajo de coverage del phase3. Si asumimos que los datos de Byrska-Bishop tienen un coverage 30X de media, no debería de ser necesario aplicarla."
		- I understand that If we have 30X coverage, we should be ok even those regions with low-depth for the low-coverage dataset, right? Or they could be still low-depth even having for 30X average depth?

- Comment jesus
	- mail 1
		- Instalación: Siguiendo tus instrucciones, he sido muy cuidadoso con el lugar de instalación del caché y los plugins, además de usar la misma versión de caché y de VEP. A la hora de instalar VEP, no me descargo el caché para evitar tener que descargarlo cada vez que construya el container. Lo que he hecho ha sido descargarlo por mi cuenta una vez (vigilando la versión) en una localización concreta que luego indico a VEP con "--dir_cache". Lo mismo he hecho con los fasta ancestrales. Mi pregunta: ¿es necesario seleccional la opción "f" (FASTA) cuando instalamos VEP con INSTALL.pl? Yo he usado "--auto ap", así que no he descargado nada más durante la instalación a parte de la api y el AncestralAllele plugin. Según el manual, "f" instala fastas que se pueden usar para incluir anotaciones HGVS, comparar con la secuencia de referencia y construir transcript models desde un GFF file. Entiendo que si yo no necesito nada de eso, no me hacen falta los fasta. De hecho, ya he corrido VEP en todos los cromosomas y no me ha dado problemas. ¿Tu recomendarías utilizar la opción "f" aun así para mi caso?
			- /home/dftortosa/singularity/dating_climate_adaptation/hg38_mig/scripts/recipes
			- http://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#fasta
		- Warnings VEP 
			- A la hora de correr INSTALL.pl, el perl script hace automaticamente varios checks para ver si la instalación se ha hecho bien. El resultado final del test en mi caso es PASS, es decir, no hay ningún error. Sin embargo, hay unos warnings que no termino de entender muy bien. Por ejemplo, cuando corre "AnnotationSource_File_GFF", dice "WARNING: The feature_type chromosome is being skipped", y lo mismo con "biological_region" , "five_prime_UTR", etc.... ¿Has visto alguna vez esto? No he encontrado información al respecto, ni hilos de github. Imagino que no será grave porque todo parece funcionar.... Adjunto el output con los checks de VEP por si le quieres echar un vistazo.
			- También obtengo este warning cada vez que corro VEP: "Smartmatch is experimental at /opt/ensembl-vep/modules/Bio/EnsEMBL/VEP/AnnotationSource/File.pm". Entiendo que esto es un warning que viene directamente de perl por usar "Smartmatch" y no de VEP. Según he leido, Smartmatch es muy problemático ya que puede dar lugar a errores, pero entiendo que no podemos hacer nada. ¿A ti te sale el mismo warning?
				- https://stackoverflow.com/questions/55819998/smartmatch-from-5-27-vs
		- Número de SNPs con/sin alelo ancestral
			- En total, hay 4,254,500 SNPs que no tienen alelo ancestral, mientras 57,344,650 si lo tienen. ¿Te parece un número razonable de missing? ¿Ó crees que mis análisis pueden estar fallando?
			- Dentro de los SNPs con alelo ancestral:
				- Hay 53,023,522 SNPs con alelo de alta confianza (apoyado por dos secuencias) y 4,321,128 con alelo de baja confianza (apoyado solo por una de las comparaciones). Esto es un 7.5%. ¿Tu usarías estos alelos con baja confianza ó directamente los eliminarías? 
				- Dentro de los SNPs con alelo ancestral (high and low confidence), hay 296,703 SNPs para los que el alelo ancestral no es el REF ni el ALT. Así que estos también serán considerados missing. ¿Te parece un número razonable ó podría ser indicativo de que hay algún problema? Imagino que no debería haber problemas con la hebra ya que Byrska-Bishop et al. hicieron strand checks comparando con el panel de referencia, así que no se a que se debe este mismatch.
		- Procesado de los alelos ancestrales: He usado bcftools +split-vep para extraer "AA" del field CSQ. Luego me he llevado "AA" a un tsv indexado con tabix donde todos los alelos ancestrales pasan a estar en mayuscula. Entonces he usado ese tsv para crear un nuevo campo con bcftools annotation. Así que al final tengo dos campos, uno con AA directamente obtenido de VEP y otro donde todos los alelos ancestrales están en mayuscula. Por tanto, el segundo campo no discrimina entre alelos de baja o alta confianza, por si al final queremos usarlos todos. ¿Le ves sentido al approach ó se te ocurre algo que me pueda faltar?
		- Máscaras de accesibilidad: A la hora de aplicar la máscara, estoy usando directamente el bed file, el cual contiene solo sitios marcados como "passed" (P). Por tanto, dejo fuera sitios marcados como L, H, Z o Q en los fastas (README). ¿Son estos los sitios que tu usas ó incluyes tambien otros marcados con LHZQ? Entiendo que al descartar LHZQ, evitamos regiones con máyor probabilidad de falsos positivos y eso es lo que queremos.
			- David say to avoid regions with a very low depth when it is not possible to do variant call.
				- Note that the average depth in past versions was 3X, but now the average is 10 times higher.
			- http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/README.accessible_genome_mask.20160622
	
	- mail 2
		- si no recuerdo mal la opción -f simplemente te descarga los fastas de referencia, pero me parece que cuando seleccionas la opción auto descarga casi todo lo necesario. Yo normalmente he seguido las instrucciones de la web y he utilizado este comando: INSTALL.pl -a cf -s homo_sapiens -y GRCh38
		- Respecto a los warnings de la instalación 
			- nunca me fijé en si me saltaban esos errores, pero puedo imaginar que como dices no tiene importancia. Eso sí, nunca me salió ningún problema relacionado con Smartmatch. En aramis, yo tengo instalado el cache de VEP (versión 109, GRCh38) en este path /labstorage/jmurgamoreno/mk_data/raw_data/homo_sapiens/vcf/vep_data/. 
			- Si no me equivoco debes de poder usarlo, por lo que quizás te recomendaría usar la imagen de singularity por defecto de VEP en el chr22 y comprobar los resultados.
				- singularity pull --name vep.sif docker://ensemblorg/ensembl-vep
				- Para poder acceder a la partición desde singularity debes de usar un punto de acceso. Lo más sencillo es usar el mismo path y ya está
				- singularity exec --bind /labstorage/jmurgamoreno/mk_data/raw_data/:/labstorage/jmurgamoreno/mk_data/raw_data/ vep.sif <comando>
		- Para extraer los ancestrales yo primero suelo hacer que haya sólo un tipo de anotación por línea y luego bcftools +split-vep. Te envío el script que uso para separar las anotaciones, pero creo que basta con bcftools.
		- Respecto a la máscara, creo que no hace falta usarla. La máscara está hecha para evitar los sitios con bajo de coverage del phase3. Si asumimos que los datos de Byrska-Bishop tienen un coverage 30X de media, no debería de ser necesario aplicarla.

	- mail 3
		- Instalación y warnings:
			- He vuelto a instalar VEP incluyendo está vez, no sólo caché, api, y plugins, sino también los fastas de referencia (opción "f"). He calculado los alelos ancestrales a lo largo de todos los cromosomas y el resultado exactamente el mismo (comparado byte by byte con cmp).
			- Respecto al warning de smartmatch, creo que ya sé qué ha pasado. El 28 de Junio de 2023 añadieron un "if" en File.pm que es lo que está causando el warning (línea 471; commit). Creo que está relacionado con variantes estructurales, porque el commit donde se añadió dice "Integrate SV Overlap plugin features". Este cambio se añadió después de publicar la versión 109 de VEP, por eso a ti no te sale pero a mí sí (yo uso la 110; historial). De hecho, para la versión 111 han propuesto cambiar esta línea con la siguiente justificación: "Smartmatch is an experimental feature that throws warnings in newer versions of Perl. This PR intends to replace a smart match with an equivalent if statement" (línea 472; link). Todo esto cuadra con el warning que me sale a mí: "Smartmatch is experimental at /opt/ensembl-vep/modules/Bio/EnsEMBL/VEP/AnnotationSource/File.pm line 472.". Por tanto, entiendo que no es un problema grave (si no lo habrían quitado de urgencia de la versión 110 nada más descubrir esto) y, en cualquier caso, sólo afectaría a variantes estructurales, no a SNPs.
				- https://github.com/Ensembl/ensembl-vep/commit/f22db24928ee96002602077b0e65751fc61fcd9b
				- https://github.com/Ensembl/ensembl-vep/commits/release/110/modules/Bio/EnsEMBL/VEP/AnnotationSource/File.pm
				- https://github.com/Ensembl/ensembl-vep/pull/1461 
			- Por todo esto, creo que mi instalación no debería ser problemática. ¿Le ves sentido? Si te sigue pareciendo raro, entonces hago la comprobación que me has dicho con la imagen por defecto en el cromosoma 22.
		- Procesado de los alelos ancestrales:
			- Siguiendo tu consejo, he hecho un awk script para directamente extraer 1 único alelo del AA field y luego intercambiar REF por ALT si el alelo ancestral es igual a ALT pero no a REF. Lógicamente, he descartado los casos donde ni REF ni ALT se corresponden con el alelo ancestral. Solo una pregunta más sobre esto, ¿tu sueles usar los alelos ancestrales con baja confianza? En los fasta con los alelos ancestrales, hay alelos en mayúscula y en minúscula, siendo low-confidence aquellos en minúscula. Se supone que en estos casos solo una de las dos comparaciones apoya que ese alelo es ancestral. David me dijo de eliminarlos si son muy pocos (1-2%), pero son 4 millones... (un 7% de todos los SNPs con alelo ancestral).
		- Máscara:
			- Ok, tiene sentido! La fase 3 tenía un coverage promedio de 7X, pero ahora tenemos 30X. De hecho, la máscara clasifica regiones según si la cobertura está por debajo o por encima de esa media de 7X. Se lo comentaré a David para asegurar que estamos todos de acuerdo, pero en principio no voy a usar la máscara.





- calculate map file using genetic distance from decode hg38 original
	- check the original map is hg38
		- it says it is hg38 in the origina file, aau1043_DataS3.gz
	- check if the data is 1 or 0 based

	- DAVID: according to DAvid, you can check that with the assembly file, if the base the map is saying is in position X is indeed at position X, you are good.



- ask david about hgdp + 1000kgp
	- A new preprint has been posted about the harmonization of 1000KP high coverage and the HGDP. 
		- A harmonized public resource of deeply sequenced diverse human genomes
		- https://www.biorxiv.org/content/10.1101/2023.01.23.525248v2.full
	- we could gain populations at higher latitudes in eurasia, which could be useful for the climate adaptation paper. 
	- although I do not know if we would need new simulations for these pops...
	- do you think it is worth it?
		- For the future could be interesting, but FIRST we need to know if flex sweep works well with 1000GP and its number of indidividuals
		- Then, in the future, we can check whether Flexsweep works with smaller population sizes, as the HGDP has less samples in many populations.
		- Months later, David sent to you again this paper! He forgot that we talked about this, and he thought that it could be interesting for me! So in the future it would be worth it to mention it and check if flex-sweep works with this smaller sample size.