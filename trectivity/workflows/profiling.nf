include { gffquant_flow } from "../../nevermore/workflows/gffquant"
include { motus3; motus4 } from "../modules/profilers/motus"


workflow profiling {

	take:
	reads_ch
	samples_ch

	main:

	if (params.run_gffquant) {
		// gq_input_ch = nevermore_main.out.fastqs
		gq_input_ch = reads_ch
			.map { sample, fastqs ->
			sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
			return [ sample_id, [fastqs].flatten() ]
		}
		.groupTuple(size: 2, remainder: true)
		.map { sample_id, fastqs -> [ sample_id, [fastqs].flatten() ] }
		.join(
			// prep_samples_ch.map { it -> [ it[0].id, it ] }, by: 0
			samples_ch.map { it -> [ it[0].id, it ] }, by: 0
		)
		// .map { sample_id, fastqs, _sample, _source, _reads, _contigs, _genes, biome -> [ sample_id, fastqs, file("${params.gq_db}/${biome}/${biome}.mmi") ] }
		// .map { sample_id, fastqs, sample_meta, _source, _reads, _contigs, _genes -> [ sample_id, fastqs, "${params.gq_db}/${sample_meta.biome}/${sample_meta.biome}.mmi" ] }
		.map { sample_id, fastqs, sample_meta -> [ sample_id, fastqs, "${params.gq_db}/${sample_meta[0].biome}/${sample_meta[0].biome}.mmi" ] }
		
		gq_input_ch.dump(pretty: true, tag: "gq_input_ch")
	
		gffquant_flow(gq_input_ch)
	}

	if (params.run_motus) {
		def run_motus3 = (params.run_motus == "motus3" || params.run_motus == "both")
		def run_motus4 = (params.run_motus == "motus4" || params.run_motus == "both")

		if (run_motus3) {
			motus3(reads_ch, params.motus3_db)
		} 
		if (run_motus4) {
			motus4(reads_ch, params.motus4_db)
		}
	}


}