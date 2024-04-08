# Things to do about hg38 mig

- CHECK ANSWER OF DAVID ABOUT MASK AND GEN POS ZERO


- final mails about hg38 mig with David
	- mail 1
		- Hi David, As you have probably seen in the email to Addison, I have already calculated the hap and map files for all the populations and chromosomes of the 1KGP. I am attaching a plot with the number of SNPs retained (upper panel) and the percentage lost (lower panel) for each chromosome and population after filtering. Chromosomes are color coded, while data points of each population have a distinct shape, showing the population names in the X axis. As expected, the number of SNPs decreased from chrom 1 to 22 and from african vs non-african populations. The final number of SNPs ranges from 120K to 1.8 million per chromosome. I guess non-african lose less SNPs just by their higher diversity, having less monomorphic SNPs. As I mentioned in my previous email, I have done many many checks in order to avoid silly errors like not selecting the correct samples or SNPs, but I remember you told me something about comparing with SNPs of the previous hg19 panel. Do you want me to do that additional check? If so, maybe, we can meet so you can give a bit of guidance about what to compare. 
	- mail 2
		- Hi Diego,  I am going to check the coordinates with a quick liftOver and let you know if everything is fine.
	- mail 3
		- Hi David, 
			- Just a few things about the new data:
			- Coordinates in the map files
        		- The map files use the selscan format, so columns are: chromosome, locus ID, genetic position, physical position.
        		- I initially calculated the map files for the 22 chromosomes with decode hg38 and then used them as templates for the map files of each population. In case it is easier for you to check the coordinates in only 22 files, you can find them here:
            		- /xdisk/denard/dftortosa/climate_adaptation_met_genes/hg38_mig/results/00b_map_files/
            		- The files are named as follows: chrXXX_selscan.map.gz
		        - Please, for future analyses, use the map files of each population. These are the final files with updated SNP after applying all the filters we discussed.
		            - /xdisk/denard/dftortosa/climate_adaptation_met_genes/hg38_mig/results/03_hap_map_files/
    		- I talked again with Jesús about the accessibility masks, and he told me that, in reality, there is no need to use them. These masks were made for the low-coverage sites of phase 3, but assuming that the new data has a 30X average depth, this is not really needed. Thus, I have not applied the masks.
		    - Genetic position. I have a few SNPs with a genetic position of 0.0 cM. 
		        - For example, in chromosome 10, the first SNP at position 616253 has 0.0 cM, while the second SNP, which is 179 bases away, has a cM value of 5.3e-43, pretty low but not zero.
		        - In chromosome 4 is the same, the first SNP has 0.0 cM, and the next one, which is 30 bases away, has 5.13e-109 cM. It is physically closer, hence its genetic position is closer to 0.
		        - This seems to be a matter of decimals (maybe decimal limit in python?) and just affecting the first SNP. Also, I visually inspected the correlation between genetic position (calculated by me) and the recombination rates (reported by decode hg38) for all chromosomes and everything seems to be fine (larger increases in genetic position under higher recombination rates). So I guess we should not have problems with this.
			- If all of this is ok for you, and you do not find anything wrong with the coordinates, I think we are good to go and the lab can start using the new hap and map files.

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
				- singularity exec --bind /labstorage/jmurgamoreno/mk_data/raw_data/:/labstorage/jmurgamoreno/mk_data/raw_data/ vep.sif
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

	- mail 4
		- disculpa la tardanza, he estado revisando y haciendo un par de scripts. Ayer además tuve un poco de fiebre.
		- Respecto al Smartmatch no me preocuparía más. Es un warning de una dependencia de VEP. Como indicas lo corregirán en la próxima versión.
		- Te envío un script de bash que instala todas las dependencias VEP desde singularity y anota los ancestrales. No obstante, he estado revisando los VCF durante el proceso y he visto que los multialélicos están en diferentes líneas. En mi caso, esto implica que al filtrar bialélicos las posiciones no se excluyen. Por ejemplo la posición 10522570 del cromosoma 22 la encontramos así:
		- chr22    10522570    22:10522570:C:A    C    A
		- chr22    10522570    22:10522571:C:T    C    T
		- Yo normalmente utilizo el siguiente comando para filtrar bialélicos.
		- bcftools view -m2 -M2 -v snps
		- Este comando no funciona por que el programa detecta que cada línea como una única mutación, no como un multialélico.
		- En el script verás una combinación de bcftools que une los multialélicos en una sola entrada y después filtra bialélicos. Todo esto te lo comento para que tengas muy en cuenta este tipo de posiciones.
		- Una vez filtrado y anotados puedes proceder a polarizar las posiciones. Yo normalmente también he incluido la annotación de baja calidad.
		- Te envío un script the python y unos comandos de plink que harían esta función. Verás que plink necesita un fichero concreto para con las posiciones y los ID para recodificar las posiciones. En mi caso verás que es un comando de awk que revisa que  el ancestral sea REF o ALT. Este comando no actualiza ningún campo del VCF, eso tendrías que hacerlo con vcftools. Desde plink puedes exportar directamente a .hap. El script de python lo hice un poco rápido para enviartelo y tampoco actualiza la INFO del VCF, revísalo bien si lo vas a usar. La única vez que yo hice esto lo hice con plink.
	- mail 4
		- Muchísimas gracias por mandarme los script. Al revisarlos he descubierto que la iba a cagar sobremanera... Yo estaba simplemente cambiando los alelos en las columnas REF y ALT, pero no en GT!! He modificado mi awk script para cambiar 1 por 0 y 0 por 1 en los genotipos de SNPs no polarizados. Además, he comprobado que obtengo los mismos resultados usando tu approach en plink2. Con esto y con la actualización de INFO usando +fill-tags parece que ya estaría solventado el problema.
		- Por cierto, aquí he detectado que aunque tu le des a plink2 una lista filtrada solo con los REF de los SNPs que tienen ancestral, luego el VCF resultante tienen todos los SNP, i.e., no filtra. También he visto que plink cambia "chr1" por "1". Yo he tenido que procesar el VCF después de plink, filtrando SNPs y añadiendo "chr" en CHROM. Te lo comento por si lo vuelves a usar y te sirve.
		- Respecto a los multialélicos, efectivamente, están separados, eso al menos no se me escapó jeje Yo estoy haciendo como dices, combinar y luego filtrar con -M y -m. También he añadido un paso anterior (probablemente poco relevante) que es eliminar SNPs monomórficos en cada población, de esa manera, si una de las filas de un SNP multialélico solo tiene REF, esa fila se perdería y el SNP pasaría a ser bialélico, por lo que se retendría al filtrar.




- calculate map file using genetic distance from decode hg38 original
	- check the original map is hg38
		- it says it is hg38 in the origina file, aau1043_DataS3.gz
	- check if the data is 1 or 0 based
		- this was answered by the first author of the paper, saying these are 1-based coordinates.
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
