include { stringtie; extract_stringtie_transcripts } from "../modules/assembly/stringtie"
include { metaT_megahit; bwa_index; bwa2assembly } from "../modules/assembly/megahit"
include { metaT_trinity } from "../modules/assembly/trinity"
include { cd_hit_est } from "../modules/assembly/cdhit"
include { quast } from "../modules/assembly/quast"


workflow assembly {
	take:
		prep_samples_ch
		alignments_ch

	main:

		stringtie(alignments_ch)

		extract_stringtie_transcripts(
			stringtie.out.gtf
				.map { sample, gtf -> [ sample.id, sample, gtf ] }
				.join(
					prep_samples_ch.map { sample -> [sample[0].id, sample[3] ] },
					by: 0
				)
				.map { sample_id, sample, gtf, genome_fasta -> [ sample, gtf, genome_fasta ] }
		)

		downstream_fq_ch = prep_samples_ch.map { meta, source, reads, contigs, genes -> [ meta, reads ] }

		assembly_input_ch = downstream_fq_ch
			.map { sample, fastqs -> 
				def meta = sample.clone()
				meta.id = meta.id.replaceAll(/\.(singles|orphans)$/, "")
				return [ meta.id, meta.sample_id, fastqs ]
			}
			.groupTuple(by: [0, 1], size: 2, remainder: true)
			.map { sample_protocol_id, sample_id, fastqs -> 
				def meta = [:]
				meta.id = sample_protocol_id
				meta.sample_id = sample_id
				return [ meta, [fastqs].flatten() ]
			}

		assembly_input_ch.dump(pretty: true, tag: "assembly_input_ch")

		metaT_megahit(assembly_input_ch, "stage1")

		metaT_trinity(assembly_input_ch, "stage1")

		quast(metaT_megahit.out.contigs.mix(metaT_trinity.out.contigs))

		cd_hit_est(
			metaT_megahit.out.contigs
				.mix(metaT_trinity.out.contigs)
				.mix(extract_stringtie_transcripts.out.transcripts)
				.map { sample, file -> [ sample.id, sample, file ] }
				.groupTuple(by: 0, size: 3)
				.map { sample_id, sample, files -> [ sample[0], files ] }
		)

		bwa_index(cd_hit_est.out.contigs)
			// metaT_megahit.out.contigs
			// 	.map { sample, contigs ->
			// 		def meta = sample.clone()
			// 		sample.assembler = "megahit"
			// 		return [ sample, contigs ]
			// 	}
			// 	.mix(
			// 		metaT_trinity.out.contigs
			// 			.map { sample, contigs ->
			// 				def meta = sample.clone()
			// 				meta.assembler = "trinity"
			// 				return [ sample, contigs ]
			// 			}
			// 	)
			// 	.mix(
			// 		extract_stringtie_transcripts.out.transcripts
			// 			.map { sample, contigs ->
			// 				def meta = sample.clone()
			// 				meta.assembler = "stringtie"
			// 				return [ sample, contigs ]
			// 			}
			// 	)

		bwa2assembly(
			downstream_fq_ch
				.map { sample, fastqs -> [ sample.id.replaceAll(/\.singles$/, ""), sample, fastqs ] }
				.combine(bwa_index.out.index, by: 0)
				.map { sample_id, sample, fastqs, index -> 
					def meta = sample.clone()
					meta.index_id = sample_id
					return [ meta, fastqs, index ]
				}
		)

		contigs_ch = cd_hit_est.out.contigs

	
	emit:
		contigs = contigs_ch
}