build_bib_file <- function() {

  # Overwrite existing file so stale entries do not persist
  bib_file <- here("FisheriesResearch_manuscript/references.bib")
  if (file.exists(bib_file)) {
    file.remove(bib_file)
  }

  # Sources with DOI
  dois <- c(
    "10.1146/annurev-statistics-031017-100427", # Aeberhard et al. 2018
    "10.1002/ecm.1470", # Auger-Methe et al. 2021
    "10.1093/icesjms/fsw046", # Berg and Nielsen, 2016
    "10.1016/j.fishres.2019.105411", # Best and Punt, 2020
    "10.18637/jss.v070.i05", # Kristensen et al., 2016
    "10.1016/j.fishres.2012.10.012", # Methot and Wetzel, 2013
    "10.1111/2041-210X.12681", # Monnahan et al., 2017
    "10.1371/journal.pone.0197954", # Monnahan and Kristensen, 2018
    "10.1093/icesjms/fsz059", # Monnahan et al., 2019
    "10.1214/088342307000000014", # Papaspiliopoulos et al., 2007
    "10.1007/978-3-031-70584-7", # Shumway and Stoffer, 2025
    #"10.1093/icesjms/fsag190", # Stawitz et al., 2027
    "10.1139/f2011-092", # Methot and Taylor, 2011
    "10.1016/j.fishres.2014.01.014", # Nielsen and Berg, 2014
    "10.1093/icesjms/fsu213", # Thorson and Minto, 2015
    "10.1139/cjfas-2013-0645", # Thorson, Jensen and Zipkin, 2014
    "10.1111/j.1467-9868.2011.00777.x", # Lindgren, Rue and Lindstrom, 2011
    "10.1016/j.fishres.2021.105967", # Stock and Miller, 2021 (WHAM)
    "10.1016/j.fishres.2021.105873", # Stock et al., 2021
    "10.1016/j.fishres.2018.10.013", # Thorson, 2019
    "10.1016/j.fishres.2024.107024", # Monnahan, 2024
    "10.1016/j.csda.2006.03.005", # Skaug and Fournier, 2006
    "10.1080/10556788.2011.597854", # Fournier et al., 2012
    "10.1139/f82-157", # Fournier and Archibald, 1982
    "10.1139/f03-104", # Maunder and Deriso, 2003
    "10.18637/jss.v076.i01", # Carpenter et al., 2017
    "10.1139/cjfas-2015-0047", # Cadigan, 2016
    "10.1139/cjfas-2017-0035", # Miller and Hyun, 2018
    "10.1093/icesjms/fsad138", # Fisch et al., 2023
    "10.1093/icesjms/fsv229", # Szuwalski and Hollowed, 2016
    "10.1073/pnas.1214879110", # Vert-pre et al., 2013
    "10.1093/icesjms/fst057", # Punt et al., 2014
    "10.1111/insr.12534", # Osgood-Zimmerman and Wakefield, 2023
    "10.1016/j.fishres.2013.09.012", # Stewart and Martell, 2014
    "10.1016/j.fishres.2013.11.001", # Martell and Stewart, 2014
    "10.1016/j.fishres.2023.106755", # Cheng et al., 2023
    "10.1139/cjfas-2017-0446", # Xu et al., 2019
    "10.1016/j.fishres.2020.105617", # Punt et al., 2020
    "10.1016/j.fishres.2012.07.025", # Maunder and Punt, 2013
    "10.1093/icesjms/fst211", # Thorson, Hicks and Methot, 2015
    "10.1139/cjfas-2015-0339", # Miller, Hare and Alade, 2016
    "10.1093/icesjms/fsu055", # Johnson et al., 2015
    "10.1093/icesjms/fsx159", # Szuwalski et al., 2018
    "10.1016/j.fishres.2013.03.015", # Deroba and Schueller, 2013
    "10.1139/cjfas-2025-0092", # Li et al., 2026
    "10.1111/faf.70082" # Cheng et al., 2026 (SPoRC)
  )
  doi_entries <- RefManageR::GetBibEntryWithDOI(dois)

  # RefManageR keys are <first author>_<year>. Several entries share a key, so
  # set those keys from the DOI rather than relying on the automatic suffix.
  fixed_keys <- c(
    "10.1093/icesjms/fsu213" = "Thorson_2014", # Thorson and Minto, 2015
    "10.1139/cjfas-2013-0645" = "Thorson_2014b", # Thorson, Jensen and Zipkin, 2014
    "10.1093/icesjms/fst211" = "Thorson_2014c", # Thorson, Hicks and Methot, 2015
    "10.1016/j.fishres.2021.105967" = "Stock_2021", # Stock and Miller, 2021
    "10.1016/j.fishres.2021.105873" = "Stock_2021b", # Stock et al., 2021
    "10.1111/faf.70082" = "Cheng_2026" # Cheng et al., 2026 (SPoRC)
  )
  entry_dois <- tolower(vapply(seq_along(doi_entries), function(i) {
    d <- doi_entries[[i]]$doi
    if (is.null(d)) NA_character_ else d
  }, character(1)))
  for (d in names(fixed_keys)) {
    idx <- which(entry_dois == tolower(d))
    if (length(idx) == 1) names(doi_entries)[idx] <- fixed_keys[[d]]
  }

  # Manual Entries
  manual_entries <- c(
    RefManageR::BibEntry(
      bibtype = "InProceedings",
      key = "Gorinova_2020",
      title = "Automatic Reparameterisation of Probabilistic Programs",
      author = "Gorinova, Maria and Moore, Dave and Hoffman, Matthew",
      booktitle = "Proceedings of the 37th International Conference on Machine Learning",
      pages = "3648--3657",
      year = "2020",
      editor = "III, Hal Daumé and Singh, Aarti",
      volume = "119",
      series = "Proceedings of Machine Learning Research",
      publisher = "PMLR",
      url = "https://proceedings.mlr.press/v119/gorinova20a.html"
    ),
    RefManageR::BibEntry(
      bibtype = "InCollection",
      key = "Betancourt_2019",
      title = "Hamiltonian Monte Carlo for hierarchical models",
      author = "Betancourt, Michael and Girolami, Mark",
      booktitle = "Current trends in Bayesian methodology with applications",
      pages = "119--142",
      year = "2019",
      publisher = "Chapman and Hall/CRC"
    ),
    RefManageR::BibEntry(
      bibtype = "Manual",
      key = "Miller_2020",
      title = "The Woods Hole Assessment Model (WHAM)",
      author = "Miller, Timothy J. and Stock, Brian C.",
      year = "2020",
      note = "Version 1.0",
      url = "https://timjmiller.github.io/wham/"
    ),
    RefManageR::BibEntry(
      bibtype = "Manual",
      key = "Hester_2025",
      title = "bench: High Precision Timing of R Expressions",
      author = "Hester, Jim and Vaughan, Davis",
      year = "2025",
      note = "R package version 1.1.4, https://github.com/r-lib/bench",
      url = "https://bench.r-lib.org/"
    ),
    RefManageR::BibEntry(
      bibtype = "Book",
      key = "Rue_2005",
      title = "Gaussian {M}arkov Random Fields: Theory and Applications",
      author = "Rue, Håvard and Held, Leonhard",
      publisher = "{Chapman and Hall/CRC}",
      address = "Boca Raton",
      year = "2005"
    ),
    RefManageR::BibEntry(
      bibtype = "Article",
      key = "Hoffman_2014",
      title = "The No-U-Turn sampler: adaptively setting path lengths in Hamiltonian Monte Carlo",
      author = "Hoffman, Matthew D. and Gelman, Andrew",
      journal = "Journal of Machine Learning Research",
      volume = "15",
      pages = "1593--1623",
      year = "2014"
    ),
    RefManageR::BibEntry(
      bibtype = "Manual",
      key = "Kristensen_2025",
      title = "RTMB: 'R' Bindings for 'TMB'",
      author = "Kristensen, Kasper",
      year = "2025",
      note = "R package version 1.9",
      url = "https://CRAN.R-project.org/package=RTMB"
    ),
    RefManageR::BibEntry(
      bibtype = "Article",
      key = "Kristensen_2026",
      title = "{RTMB}: An {R} Package for Fitting Latent Variable Models With {T}emplate {M}odel {B}uilder",
      author = "Kristensen, Kasper and Nielsen, Anders and Albertsen, Christoffer Moesgaard and Berg, Casper Willestofte and Trijoulet, Vanessa",
      journal = "Ecology and Evolution",
      volume = "16",
      number = "11",
      pages = "e74572",
      year = "2026",
      doi = "10.1002/ece3.74572"
    ),
    RefManageR::BibEntry(
      bibtype = "Book",
      key = "Quinn_1999",
      title = "Quantitative Fish Dynamics",
      author = "Quinn, Terrance J. and Deriso, Richard B.",
      publisher = "Oxford University Press",
      address = "New York",
      year = "1999"
    ),
    RefManageR::BibEntry(
      bibtype = "Misc",
      key = "Monnahan_2026",
      title = "Leveraging sparsity to improve no-{U}-turn sampling efficiency for hierarchical {B}ayesian models",
      author = "Monnahan, Cole C. and Kristensen, Kasper and Thorson, James T. and Carpenter, Bob",
      howpublished = "arXiv preprint arXiv:2603.02437",
      url = "https://arxiv.org/abs/2603.02437",
      year = "2026"
    )
  )

  RefManageR::WriteBib(c(doi_entries, manual_entries),
                       file = bib_file)

  # Tidy up the written file ... some of them the pandoc cant read
  bib_txt <- readLines(bib_file, warn = FALSE)
  bib_txt <- gsub("\\a'", "\\'", bib_txt, fixed = TRUE)
  bib_txt <- gsub("&amp;#8217;", "'", bib_txt, fixed = TRUE)
  bib_txt <- gsub("<scp>", "{", bib_txt, fixed = TRUE)
  bib_txt <- gsub("</scp>", "}", bib_txt, fixed = TRUE)
  bib_txt <- gsub("<b>", "{", bib_txt, fixed = TRUE)
  bib_txt <- gsub("</b>", "}", bib_txt, fixed = TRUE)
  bib_txt <- gsub("<i>", "{", bib_txt, fixed = TRUE)
  bib_txt <- gsub("</i>", "}", bib_txt, fixed = TRUE)
  bib_txt <- gsub("&amp;", "\\&", bib_txt, fixed = TRUE)
  bib_txt <- gsub("([[:alpha:]}])\\?([[:alpha:]{])", "\\1-\\2", bib_txt)
  bib_txt <- paste(bib_txt, collapse = "\n")
  bib_txt <- gsub("}\n\\s*:", "}:", bib_txt)
  writeLines(bib_txt, bib_file)
}
