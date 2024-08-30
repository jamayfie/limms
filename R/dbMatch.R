#' Matches unknown ions to a database of known masses
#'
#' This function allows an XCMS diffreport object, limma topTable object, or 
#' other table to be annotated by matching the mass to an object with known 
#' masses. Retention time can also be matched to known times.
#' 
#' This function matches metabolite peaks to the database using a specified
#' mass window (+/- parts per million) and time window (+/- seconds); hence,
#' matches are to compounds within the specified window and are putative.
#' One peak may match multiple compounds in the database within the windows,
#' or multiple peaks may match the same compound in the database.  This can
#' reflect distinct compounds close in mass, compounds with the same chemical
#' formula but a different atomic arrangement, and/or technical artifacts.
#' 
#' dbMatch is inclusive, meaning multiple matches are allowed and the output
#' may contain new rows for each of the matches.  Users must resolve such
#' conflicts with other/additional information.  Including chromatographic 
#' retention times can increase accuracy.  dbMatch can use retention time
#' matching in addition to mass matching, but does not require time matches.
#' The database may contain retention time data for only a subset of
#' of entries.  Unlike mass matches, which are required to generate an
#' annotation, retention time matches are coded 'Y'/'N' for yes/no.
#' 
#' Many annotations will be incorrect!  dbMatch is a discovery tool intended
#' to filter noise, guide and prioritize hits.
#' 
#' In no case are these annotations publishable without further validation.
#'
#' @param	x	An R object containing the data to be annotated.  At minimum, it
#' must include a column of masses, for example, the XCMS diffreport mzmed.
#' @param	peaknames  The name of the column in x containing the peaknames.
#' @param	mzmedcol	The name of the column in x containing masses.  This is 
#' typically the median mz from a group of measurements, which will be matched
#' to the database by extrapolating a ppm window around the known masses.
#' For an alternative approach, see the 'mzmin'/'mzmax' arguments.
#' @param	mzmin	Specify an exact minimum and maxinum mz window instead of
#' using a ppm window around the median.  Set 'mzmin' to the name of the column
#' containing the minimum mass to use this feature.
#' If 'mzmin' and 'mzmax' are not both set, a mzmedcol is required and is used.
#' If using 'xcms', 'mzmin' and 'mzmax' columns are standard output, but may vary
#' quite a bit in their ppm window size(s).
#' There may be some advantage to this approach for odd peak shapes.
#' @param	mzmax	Specify an exact minimum and maxinum mz window instead of
#' using a ppm window around the median.  Set 'mzmax' to the name of the column
#' containing the maximum mass to use this feature.  
#' @param	db	An R object containing known/named compounds and their masses.
#' Additional columns of masses for alternative adducts, 
#' retention time data, annotations, metabolite classes, etc are useable.
#' @param	dbcol	The column(s) in db containing masses.  Column names or number
#' are accepted.  Multiple columns can be specified with c(), although a single
#' column of masses with annotations for features like adducts and isotopes is 
#' preferred. If multiple columns are specified, they are collated into a
#' "db_columns" column in the output.
#' @param	ppmCut	The mass tolerance for matching in parts per million (ppm).
#' The default setting, 10 ppm, selects a window of 10 ppm below and above
#' the known mass.
#' @param	rtmedcol	The column in x containing the median retention time. Leave
#' blank if retention time matching is not desired.  As for mass, an alternative
#' to using a time window is to directly specify the minimum and maximum
#' retention time for each peak using the 'mzmin' and 'mzmax' arguments.
#' @param	rtmin	Specify an exact minimum and maxinum retention time instead of
#' using a time window around the median.  Specify the name of the column in x
#' containing the minimum retention time.
#' @param	rtmax	Specify an exact minimum and maxinum retention time instead of
#' using a time window around the median.  Specify the name of the column in x
#' containing the maximum retention time.
#' If 'rtmin' and 'rtmax' are not both set, a rtmedcol is required and is used.
#' If using 'xcms', 'rtmin' and 'rtmax' columns are standard output.
#' As for 'mzmin'/'mzmax', there may be some advantage to using the measured
#' retention time window with the caveat that the window sizes vary.
#' @param	x_rt	Can be used to allow retention times in seconds or minutes.
#' The default x_rt=60 assumes rtmedcol is in seconds and dbrt is in minutes,
#' the case for an xcms output matched to a typical database.  Set x_rt=1 if
#' the input data and database are both in seconds or both in minutes.
#' Set x_rt=(1/60) if the input data is in minutes and the database in seconds.
#' @param	dbrt	The name of the column in db containing retention times.
#' Leave blank if retention time matching is not desired.
#' @param	rt_tolerance	The desired retention time tolerance in seconds. The
#' default value is 180 seconds above and below the known retention time in db.
#' @param	db_annotations	The names of additional column(s) to append to database 
#' hits, such as chemical formulas, abbreviations, full compound names, or notes.
#' These columns are passed to the final object, but are also used in a call to
#' 'unique' to collapse redundant annotations.  Hence, to separately annotate
#' compounds with the same mass, retention time, and chemical formula, but 
#' different atomic arrangements, specify a database column that makes these
#' annotations unique (ie and ID column with rows for Leucine and Isoleucine).
#' To pass multiple columns, set db_annotations=c("column1", "column2"), etc.
#' @note dbMatch identifies windows consistent with compounds, not compounds.
#' In many cases a peak will match more than one known.  
#' This function is provided to suggest additional experiments only.
#' @details	This function provides a method to append putative names to
#' compounds with a similar mass based on matches to a user provided list
#' of known compounds and masses.  Users must provide a list of knowns
#' with masses, which may be limited to spike-in controls and calibration
#' compounds, or may be as extensive as a list of all metobolites known to be 
#' found in a certain cell type.
#'
#' Multiple columns of known masses can be specified by dbcol to allow
#' matching to adducts, or the adduct annotation can be provided in a
#' separate column.
#'
#' The db_annotations argument allows columns of additional annotation
#' to be included; for example, chemical formulas or compound classes.  
#' 
#' The dbMatch output is a data.table with annotations columns appended to
#' the input object.  The appended columns include the known mass ('exact_mz'), 
#' if using retention time windows, a column of 'Y'/'N' indicating whether 
#' or not the mass match was a retention time match ('rt_match'), a column of
#' 'adduct' data, and any additional columns passed through 'db_annotation'.
#' 
#' @import data.table
#' @importFrom gdata drop.levels
#' @export
#' @examples
#' 
#' require(data.table)
#' 
#' # Use the example CBS dataset CBS.xcms_diffreport, an xcms diffreport
#' names(CBS.xcms_diffreport)
#'
#' # Read in a list of knowns. hmdb_Shortlist contains a list
#' # of compounds used as calibration standards in the CBS experiment.
#' hmdb <- hmdb_Shortlist
#'
#' # Run dbMatch.  The mzmed column in the object all is column 6.
#' x.ID_db <- dbMatch(x=CBS.xcms_diffreport, peaknames="name", mzmedcol="mzmed", 
#'   db=hmdb, dbcol="Adduc_MW_.Da.", ppmCut=10, db_annotations=c("ID", 
#'   "Formula", "HMDB_ID", "Adduct"))
#'
#' # Both the peaklist and database have an "ID" column, so the merge failed with an error message.
#'
#' # Change the name of the database column and re-run
#' names(hmdb)[1] <- "Metabolite"
#'
#' x.ID_db <- dbMatch(x=CBS.xcms_diffreport, peaknames="name", mzmedcol="mzmed", db=hmdb, 
#'   dbcol="Adduc_MW_.Da.", ppmCut=10, 
#'   db_annotations=c("Metabolite", "Formula", "HMDB_ID", "Adduct"))
#' # A fast check: pull only the database hits
#' x.ID_db[!is.na(exact_mz), c(1, 6, 70:75)]
#'
#'
#' # dbMatch can be used on other input types, like a limma topTable object
#' # so that only significant hits are annotated.
#' # Read in a topTable, such as significant hits from the CBS dataset
#' # ranked by F score in the dataset ttF.CBS, 
#' # then merge with the mass info, keeping the topTable order
#' 
#' ttF_anno <- merge(data.table(ttF.CBS), CBS.xcms_diffreport, by.x="ID", by.y="name", all.x=TRUE, 
#'   all.y=FALSE, sort=FALSE)
#'
#' x.ID_db <- dbMatch(x=ttF_anno, peaknames="ID", mzmedcol="mzmed", db=hmdb, 
#'   dbcol=4, ppmCut=10, 
#'   db_annotations=c("Metabolite", "Formula", "HMDB_ID", "Adduct"))
#' # Now with a proper adjusted p value
#' x.ID_db[!is.na(exact_mz), c(1, 9, 14, 78:83)]


dbMatch <-
  function(x,
           peaknames = "rn",
           mzmedcol = "mz",
           mzmin = NULL,
           mzmax = NULL,
           db,
           dbcol,
           ppmCut = 10,
           rtmedcol = "rt",
           rtmin = NULL,
           rtmax = NULL,
           x_rt = 60,
           dbrt,
           rt_tolerance = 180,
           db_annotations = NULL) {
    exact_mz <- ppm.start <- ppm.end <- rt_match <- maxRt <- NULL
    rt.start <- rt.end <- mzId <- minMz <- maxMz <- . <- NULL

    if (is.data.table(x) == FALSE) {
      x <- data.table(x)
      x.peaks <- x[, ..peaknames]
    }
    if (is.data.table(db) == FALSE) {
      db <- data.table(db)
    }
    
    if (any(names(x) %in% names(db))) {
      print("Error: x and db have column(s) with the same name. Please change duplicated names.")
    }
    
    if ((length(dbcol) > 1) == TRUE) {
      db.melt <-
        data.table::melt(
          db,
          measure.vars = dbcol,
          variable.name = "db_columns",
          value.name = "exact_mz",
          na.rm = TRUE
        )
    }
    else {
      db.melt <-
        data.table::melt(
          db,
          measure.vars = dbcol,
          value.name = "exact_mz",
          na.rm = TRUE
        )
    }
    
    ppmT <- ppmCut / 1000000
    
    if (missing(dbrt)) {
      # no retention time matching
      print("No retention time detected.  Mass match only")
      
      if (is.null(mzmin) | is.null(mzmax)) {
        # using mzmed
        x[, c("minMz", "maxMz") := x[, .SD, .SDcols = mzmedcol]]
        db.range <-
          db.melt[, c("ppm.start", "ppm.end") := .(-(ppmT * exact_mz) + exact_mz, (ppmT *
                                                                                     exact_mz) + exact_mz)]
        data.table::setkey(db.range, ppm.start, ppm.end)
        
        x.db <-
          data.table::foverlaps(x,
                    db.range,
                    by.x = c("minMz", "maxMz"),
                    type = "within")
        x.db[, rt_match := "N"]
      }
      
      else {
        # using mzmin-mzmax
        x[, c("minMz", "maxMz") := .(mzmin, mzmax)]
        db.range <-
          db.melt[, c("ppm.start", "ppm.end") := .(exact_mz - 1e-5, exact_mz + 1e-5)] # 1e-5 to jitter so min != max
        data.table::setkey(db.range, ppm.start, ppm.end)
        
        x.db <-
          data.table::foverlaps(x,
                    db.range,
                    by.x = c("minMz", "maxMz"),
                    type = "within")
        x.db[, rt_match := "N"]
      }
    }
    
    else {
      print(
        "The dbMatch default is to assume 'rtmedcol' is seconds (the xcms default) and 'dbrt' is in minutes (the MycoMassDB default). dbMatch sets dbrt to seconds using 'dbrt*x_rt'; hence, correct these assumptions by resetting x_rt accordingly. Likewise, rt_tolerance uses seconds."
      )
      
      if (is.null(rtmin) | is.null(rtmax)) {
        x[, c("minRt", "maxRt") := x[, .SD, .SDcols = rtmedcol]]
        x[, maxRt := maxRt + 0.1] # min != max
      }
      
      else {
        x[, c("minRt", "maxRt") := .(rtmin, rtmax)]
        x[, maxRt := maxRt + 0.1] # min != max
      }
     
      db.melt[, c("rt.start", "rt.end") := db[, .SD, .SDcols = ..dbrt]]
      data.table::set(db.melt, which(is.na(db.melt$rt.start)), "rt.start",-1001)
      data.table::set(db.melt, which(is.na(db.melt$rt.end)), "rt.end",-1000)
      
      if (identical(db.melt$rt.start, db.melt$rt.end)) {
        #if start & end rt
        db.melt <- stats::na.omit(db.melt, cols = "rt.start")
        db.melt[, rt.start := rt.start * x_rt - rt_tolerance]
        db.melt[, rt.end := rt.end * x_rt + rt_tolerance]
      }
      
      else {
        db.melt <- db.melt[!(is.na(db.melt$rt.start) & is.na(db.melt$rt.end))] #remove NA rows
        # replace single NAs
        db.melt$rt.start[is.na(db.melt$rt.start)] <-
          db.melt$rt.end[is.na(db.melt$rt.start)]
        db.melt$rt.end[is.na(db.melt$rt.end)] <- db.melt$rt.start[is.na(db.melt$rt.end)]
        
        db.melt[, rt.start := rt.start * x_rt - rt_tolerance]
        db.melt[, rt.end := rt.end * x_rt + rt_tolerance]
      }

      
      if (is.null(mzmin) | is.null(mzmax)) {
        x[, c("minMz", "maxMz") := x[, .SD, .SDcols = mzmedcol]]
        db.range <-
          db.melt[, c("ppm.start", "ppm.end") := .(-(ppmT * exact_mz) + exact_mz, (ppmT *
                                                                                     exact_mz) + exact_mz)]
        data.table::setkey(db.range, ppm.start, ppm.end)

        x.rt <-
          data.table::foverlaps(
            x,
            db.range,
            by.x = c("minMz", "maxMz"),
            type = "within",
            nomatch = 0
          )
      }
      
      else {
        x[, c("minMz", "maxMz") := .(mzmin, mzmax)]
        db.range <-
          db.melt[, c("ppm.start", "ppm.end") := .(exact_mz - 1e-5, exact_mz + 1e-5)] # 1e-5 to jitter so min != max
        data.table::setkey(db.range, ppm.start, ppm.end)

        x.rt <-
          data.table::foverlaps(
            x,
            db.range,
            by.x = c("minMz", "maxMz"),
            type = "within",
            nomatch = 0
          )
      }

      x.rt[, mzId := paste(peaknames, exact_mz, minMz, maxMz, sep = "_")]

      rt.range <-
        x.rt[, .(mzId = mzId,
                 rt.start = rt.start,
                 rt.end = rt.end)]
      data.table::setkey(rt.range, mzId, rt.start, rt.end)

      x.db <-
        data.table::foverlaps(
          x.rt,
          rt.range,
          by.x = c("mzId", "minRt", "maxRt"),
          type = "within",
          nomatch = NA
        )
      x.db[, rt_match := ifelse(!is.na(rt.start), "Y", "N")]
    }
    
    
    x.db <- stats::na.omit(x.db, "exact_mz")
    
    if ((length(dbcol) > 1) == TRUE) {
      dbHits <-
        x.db[, .SD, .SDcols = c(peaknames,
                                "exact_mz",
                                "rt_match",
                                "db_columns",
                                db_annotations)]
    }
    else {
      dbHits <-
        x.db[, .SD, .SDcols = c(peaknames, "exact_mz", "rt_match", db_annotations)]
    }
    
    x.db <-
      merge(
        x,
        dbHits,
        by.x = peaknames,
        by.y = peaknames,
        all = TRUE,
        no.dups = FALSE
      )
    x.db <- unique(x.db)
    x.db <- gdata::drop.levels(x.db)
    x.db[, c("minMz", "maxMz") := NULL]
    if (!missing(dbrt)) {
      x.db[, c("minRt", "maxRt") := NULL]
    }
    x.db
  }
