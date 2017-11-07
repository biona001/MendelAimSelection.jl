using MendelAimSelection
using MendelBase
using SnpArrays
using Compat
import Compat: view
using DataFrames 
using Distributions

@testset "construct_pvalue_vec!" begin
    #
    # Given a SNP, this function loops through every person and keeps track of 2 vectors 
    # alleles and genes (explained below), each counting some number for a specific race. 
    # Then using these two vectors 
    # 
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    process_keywords!(keyword, "AIM 1000genomes_chr1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    populations = person.populations
    people = person.people
    snps = snpdata.snps
    ethnic = blanks(people)
    copy!(ethnic, pedigree_frame[:Ethnic])
    population = unique(ethnic)
    populations = length(population)
    alleles = zeros(populations)
    # alleles = vector of races. Given a SNP, it counts how many people in each race contains
    # that SNP in your population sample. Half will be added for those on sex chromosomes.
    genes = zeros(populations)
    # genes = vector of races. Given a SNP, it counts how many copies of the chromosome
    # under consideration is in a race (1 on each chromosome for non-sex linked genes)
    dosage = zeros(people) # a particular SNP, 1 if person has that SNP, 0 otherwise
    pvalue = ones(snps)
    
    for i in 1:5 # 1st to 5th SNP
        MendelAimSelection.construct_pvalue_vec!(dosage, snpdata, i, alleles, genes, ethnic, 
        population, person, populations, people, pvalue)
    end

    for i in 100504:100508 # last 5 SNPs
        MendelAimSelection.construct_pvalue_vec!(dosage, snpdata, i, alleles, genes, ethnic, 
        population, person, populations, people, pvalue)
    end

    @test signif(pvalue[1], 5) == 0.0021199
    @test signif(pvalue[2], 5) == 0.0022310 # would be 1 in real run because minor allele freq < 1
    @test signif(pvalue[3], 5) == 4.95e-8
    @test signif(pvalue[4], 5) == 1.7162e-5
    @test signif(pvalue[5], 5) == 0.033521

    @test signif(pvalue[100504], 5) == 0.0020991
    @test signif(pvalue[100505], 5) == 0.00082908
    @test signif(pvalue[100506], 5) == 0.00082908
    @test signif(pvalue[100507], 5) == 0.024556
    @test signif(pvalue[100508], 5) == 0.061064
end

@testset "basics" begin
    #
    # First test run aim_selection_option on 2 sample datasets, then test run 
    # the wrapper function AimSelection. aim_selection_option takes the pvalue vector
    # and ranks them. 
    #
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    process_keywords!(keyword, "AIM 1000genomes_chr1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    execution_error = MendelAimSelection.aim_selection_option(person, snpdata,
    pedigree_frame, snp_definition_frame, keyword)

    @test execution_error == false

    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    process_keywords!(keyword, "AIM 1000genomes_chrX Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    execution_error = MendelAimSelection.aim_selection_option(person, snpdata,
    pedigree_frame, snp_definition_frame, keyword)

    @test execution_error == false

    final_test_1 = AimSelection("AIM 1000genomes_chr1 Control.txt")
    final_test_2 = AimSelection("AIM 1000genomes_chrX Control.txt")

    @test final_test_1 == nothing #i.e. no error, so everything went through
    @test final_test_2 == nothing #i.e. no error, so everything went through
end





