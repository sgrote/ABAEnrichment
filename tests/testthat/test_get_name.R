

# context("get sampled substructures")

test_that("names from developmental ontology are returned",{
  	expect_that(get_name("Allen:13324"), is_equivalent_to("VLTC_ventrolateral temporal neocortex"))
  	expect_that(get_name(c("Allen:10208","Allen:10209")), is_equivalent_to(c("PCx_parietal neocortex","S1C_primary somatosensory cortex (area S1, areas 3,1,2)")))
  	expect_that(get_name(c("Allen:10209","Allen:10208")), is_equivalent_to(c("S1C_primary somatosensory cortex (area S1, areas 3,1,2)","PCx_parietal neocortex")))
	expect_that(names(get_name(c("Allen:10208","Allen:10209"))), equals(c("Allen:10208","Allen:10209")))
})

test_that("names from adult ontology are returned",{
  	expect_that(get_name("Allen:4010"), is_equivalent_to("PrG_precentral gyrus"))
  	expect_that(get_name(c("Allen:9222","Allen:9223")), is_equivalent_to(c("cc_corpus callosum", "gcc_genu of the corpus callosum")))
  	expect_that(get_name(c("Allen:9223","Allen:9222")), is_equivalent_to(c("gcc_genu of the corpus callosum","cc_corpus callosum")))
	expect_that(names(get_name(c("Allen:9222","Allen:9223"))), equals(c("Allen:9222","Allen:9223")))
})

test_that("names from different ontologies can be requested together",{
  	expect_that(get_name(c("Allen:10208","Allen:9222")), is_equivalent_to(c("PCx_parietal neocortex","cc_corpus callosum")))

})

test_that("error message for invalid structure is returned",{
	expect_that(get_name("Allen:123"), throws_error("Invalid structure_id: Allen:123."))
	expect_that(get_name(c("Allen:9223","Allen:123")), throws_error("Invalid structure_id: Allen:123."))
})


