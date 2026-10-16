params.velvet_hash_length = "25,97,4"
params.velveth_kmin = 25
params.velveth_kmax = 97
params.velveth_kstep = 4



process metaT_velvetoptimiser {
	container "quay.io/biocontainers/perl-velvetoptimiser:2.2.5--0"
	label "velvet"

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
	
	def outdir = "assemblies/metaT_velvet/${sample.id}/${stage}/${sample.library_source}/${sample.id}/${sample.id}_velveth"

	"""
	VelvetOptimiser.pl \
		-t 2 \
		-v \
		-f '${input_files}' \
		-s ${params.velveth_kmin} \
		-e ${params.velveth_kmax} \
		-x ${params.velveth_kstep} \
		-d ${outdir} \
		-p ${sample.id} \
		-o '-ins_length ${insert_size} -read_trkg yes' \
		-a

	"""


}


process metaT_velvet {
	container "quay.io/biocontainers/velvet:1.2.10--h577a1d6_9"
	label "velvet"

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