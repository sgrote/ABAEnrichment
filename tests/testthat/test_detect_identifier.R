  
# context("detect_identifier")

test_that("Entrez is recognized",{
  	expect_that(detect_identifier("34656"), equals("entrezgene"))
	expect_that(detect_identifier(34656), equals("entrezgene"))
	expect_that(detect_identifier(as.factor("34656")), equals("entrezgene"))
	expect_that(detect_identifier(as.factor(34656)), equals("entrezgene"))
})

test_that("ENSEMBL is recognized",{
  	expect_that(detect_identifier("ENSG34656918234"), equals("ensembl_gene_id"))
	expect_that(detect_identifier(as.factor("ENSG34656918234")), equals("ensembl_gene_id"))
})

test_that("ENSEMBL is not recognized",{
  	expect_that(detect_identifier("ENSG346569182")=="ensembl_gene_id", is_false())
  	expect_that(detect_identifier(as.factor("ENSG346569182"))=="ensembl_gene_id", is_false())
  	expect_that(detect_identifier("ENS346569182345")=="ensembl_gene_id", is_false())

})

test_that("HGNC is recognized",{
  	expect_that(detect_identifier("A4GNT"), equals("hgnc_symbol"))
	expect_that(detect_identifier(as.factor("A4GNT")), equals("hgnc_symbol"))
	expect_that(detect_identifier(as.factor("4AGNT")), equals("hgnc_symbol"))
	expect_that(detect_identifier(as.factor("A4GNTAH56THFGBH")), equals("hgnc_symbol"))
})
