process runGraphST {
    // Define unique tag for logs
    tag "graphst_${tag_suffix}"
    
    // Labels determine resource allocation (CPU/GPU/Time)
    label "retry" 
    label "longer_time"
    label ( params.gpu ? "use_gpu" : "use_cpu" ) // Request GPU if enabled in config

    // Point to your custom container
    container 'xinruoz/graphst_custom'
    
    // Print output to log
    echo true

    input:
        // Input: Spatial data (H5AD + RDS tuple) and Single-cell data (H5AD)
        tuple path(sp_input), path(sp_input_rds)
        path(sc_input)

    output:
        // Output: A tuple expected by the benchmark pipeline
        // (method_name, output_file, original_rds_file)
        tuple val('graphst_custom'), path("$output"), path(sp_input_rds)

    script:
        // Define file names and parameters
        tag_suffix = file(sp_input).getSimpleName()
        output = "proportions_graphst_${tag_suffix}.tsv"
        
        // Check if GPU is requested in nextflow config
        cuda_device = ( params.gpu ? "cuda" : "cpu" )
        
        // Define which model variant to use (GAT, SGC, GATv2, or 10X/default)
        // You can control this via --model_type in the Nextflow command line if you add it to params
        model_type = params.graphst_model_type ?: '10X' 

        """
        # Command to run your python script
        python $baseDir/subworkflows/deconvolution/GraphST/script_nf.py \\
            $sc_input \\
            $sp_input \\
            $params.annot \\
            --datatype $model_type \\
            --device $cuda_device \\
            --output $output
        """
}