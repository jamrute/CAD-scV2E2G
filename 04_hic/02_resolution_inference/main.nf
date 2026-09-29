nextflow.enable.dsl=2

workflow {

    // Read samplesheet into a channel
    samples = Channel.fromPath(params.samplesheet)
                 .splitCsv(header: true, sep: ',')
                 .map { row ->
                     def entries = []
                     assert row.containsKey('sample') : "CSV file does not have 'sample' column"
                     row.findAll { it.key.startsWith('bam') && it.value != '' }.each { key, value ->
                         def entry = [row.sample, file(value)]
                         entries.add(entry)
                     }
                     return entries
                 }
                 .flatMap { it } // Flatten the list
                 .groupTuple(by: 0) // Group by sample

    process mergeBamFiles {
        container '<AWS_account_#>.dkr.ecr.<AWS_region>.amazonaws.com/quay/biocontainers/samtools:1.17--h00cdaf9_0'

        input:
        tuple val(sample_id), path(bam_files)

        output:
        tuple val(sample_id), path("${sample_id}_merged.bam")

        cpus 48
        memory { 96.GB }

        script:
        """
        samtools merge -@ ${task.cpus} ${sample_id}_merged.bam ${bam_files.join(' ')}
        """
    }

    process filterBamFile {
        container '<AWS_account_#>.dkr.ecr.<AWS_region>.amazonaws.com/hicup:0.9.2--hdfd78af_1'

        input:
        tuple val(sample_id), path(bam_file)
        path genome_digest

        output:
        tuple val(sample_id), path("${bam_file.baseName}.filt.bam"), emit: bam
        path("${sample_id}_QC_outputs/*"), emit: qc_outputs

        cpus 48
        memory { 96.GB }

        publishDir '<path>/results/QC_outputs', pattern: "${sample_id}_QC_outputs/*", mode: 'copy'

        script:
        """
        hicup_filter --digest ${genome_digest} --longest 800 --shortest 150 --zip --threads ${task.cpus} ${bam_file}
        mkdir ${sample_id}_QC_outputs
        mv hicup_filter_* ${sample_id}_QC_outputs/
        """
    }

    process runHiCRes {
        container '<AWS_account_#>.dkr.ecr.<AWS_region>.amazonaws.com/hicres:cli'

        input:
        tuple val(sample_id), path(bam_file)
        path chromsizes_file

        output:
        path("${sample_id}_HiCRes_outputs/*")

        cpus 96
        memory { 192.GB }

        publishDir '<path>/results/', mode: 'copy'

        script:
        """
        mv ${bam_file} /tmp/
        mv ${chromsizes_file} /tmp/
        hicres -m bam -c ${chromsizes_file.name} -t 96 -b ${bam_file.name}
        mkdir ${sample_id}_HiCRes_outputs
        cd ${sample_id}_HiCRes_outputs
        mv /tmp/hicres/* ./
        """
    }

    //Create channels for genome digest and chromsizes_file
    ch_genome_digest = file(params.genome_digest)
    ch_chromsizes_file = file(params.chromsizes_file)

    //Run the pipeline!
    mergeBamFiles(samples)
    filterBamFile(mergeBamFiles.out, ch_genome_digest)
    runHiCRes(filterBamFile.out.bam, ch_chromsizes_file)
}
