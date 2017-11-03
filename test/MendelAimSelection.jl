using MendelAimSelection
using MendelBase
using SnpArrays
using Compat
import Compat: view
using DataFrames 
using Distributions

@testset "basics" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    process_keywords!(keyword, "AIM 1000genomes_chr1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    
end