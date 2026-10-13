params.kmer_steps = "25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93,97"

process metaT_rnaspades {
	container "quay.io/biocontainers/spades:4.3.0--hde4eca7_2"
	label "spades"

	input:
	tuple val(sample), path(fastqs)
	val(stage)

	output:
	tuple val(sample), path("assemblies/metaT_rnaspades/${stage}/${sample.library_source}/${sample.id}/${sample.id}.${stage}.*.fasta"), emit: contigs

	script:
	def mem_gb = task.memory.toGiga()
	def mem = task.memory.toBytes()

	def input_files = ""
	// we cannot auto-detect SE vs. PE-orphan!
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	if (r1_files.size() != 0) {
		input_files += "--pe1-1 ${r1_files.join(' ')}"
	}
	if (r2_files.size() != 0) {
		input_files += " --pe1-2 ${r2_files.join(' ')}"
	}
	if (orphans.size() != 0) {
		input_files += " --s1 ${orphans.join(' ')}"
	}

	script:
	def outdir = "assemblies/metaT_rnaspades/${stage}/${sample.library_source}/${sample.id}"

	"""
	rnaspades.py -t ${task.cpus} -k ${params.kmer_steps} -o ${outdir} ${input_files}
	"""



}