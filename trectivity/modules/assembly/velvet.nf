params.velvet_hash_length = "25,97,4"

process metaT_velvet {
	container "docker://quay.io/biocontainers/velvet:1.2.10--h577a1d6_9"
	label "megahit"

	input:
	tuple val(sample), path(fastqs), val(insert_size)
	val(stage)

	output:
	// tuple val(sample), path("assemblies/metaT_velvet/${stage}/${sample.library_source}/${sample.id}/${sample.id}.${stage}.*.fasta"), emit: contigs
	tuple val(sample), path("assemblies/metaT_velvet/${sample}/${stage}/${sample.library_source}/${sample.id}/${sample.id}_velveth"), emit: velveth
	
	script:
	def input_files = ""
	// we cannot auto-detect SE vs. PE-orphan!
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	if (r1_files.size() != 0 && r2_files.size() != 0) {
		input_files += "-shortPaired -fmtAuto -separate ${r1_files[0]} ${r2_files[0]}"
		if (orphans.size() != 0) {
			input_files += " -short2 -fmtAuto ${orphans[0]}"
		}
	} else {
		if (r1_files.size() != 0) {
			input_files += "-short -fmtAuto ${r1_files[0]}"
		} else if (r2_files.size() != 0) {
			input_files += "-short -fmtAuto ${r2_files[0]}"
		}
		if (orphans.size() != 0) {
			input_files += " -short2 -fmtAuto ${orphans[0]}"
		}
	}
	
	def outdir = "assemblies/metaT_velvet/${sample}/${stage}/${sample.library_source}/${sample.id}/${sample.id}_velveth"

	"""
	velveth ${outdir} ${params.velvet_hash_length} ${input_files}

	velvetg ${outdir} -ins_length ${insert_size} -read_trkg yes
	"""

}