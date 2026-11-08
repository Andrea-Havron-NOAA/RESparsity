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
    "10.1007/978-3-031-70584-7"#, # Shumway and Stoffer, 2025
    #"10.1093/icesjms/fsag190" # Stawitz et al., 2027
  )
  doi_entries <- RefManageR::GetBibEntryWithDOI(dois)

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
    )
  )

  RefManageR::WriteBib(c(doi_entries, manual_entries),
                       file = bib_file)
}
