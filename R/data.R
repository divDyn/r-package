#' Fossil occurrences of scleractinian (stony) corals from the Paleobiology Database
#' 
#' Example dataset to illustrate the package's basic functionalities.
#'
#' This particular dataset was used in a study by Kiessling and Kocsis (2015). All occurrences of Scleractinia were downloaded from the Paleobiology Database (PaleoDB, \url{http://www.paleobiodb.org/}) on 23 September 2014, originally comprising 32420 occurrences. They were than cross-checked with data from Corallosphere (\url{http://corallosphere.org}). See the article text for details.
#' 
#' References
#' 
#' Kiessling, W., & Aberhan, M. (2007). Environmental determinants of marine benthic biodiversity dynamics through Triassic–Jurassic time. Paleobiology, 33(3), 414-434.
#' 
#' @format A \code{data.frame} with 29775 observations and 38 variables:
#' \describe{
#' 	\item{\code{genus}}{Genus names of the occurrences. Cross referenced with a compiled table, the simplified version of this can be found in the supplementary material of Kiessling and Kocsis (2015).}
#' 	\item{\code{collection_no}}{The number of the collection of the occurrence in the PaleoDB.}
#' 	\item{\code{family}}{Family name of the occurrence.}
#' 	\item{\code{abund_value}}{Abundance value.}
#' 	\item{\code{abund_unit}}{Unit of abundance values.}
#' 	\item{\code{reference_no}}{The reference number of the occurrence in the PaleoDB.}
#' 	\item{\code{life_habit}}{The lifestyle of the occurring taxon.}
#' 	\item{\code{diet}}{The diet of the occurring taxon.}
#' 	\item{\code{country}}{Country of occurrence.}
#' 	\item{\code{geoplate}}{Plate id of the occurrence.}
#' 	\item{\code{lat}}{Present day latitude of the occurrence.}
#' 	\item{\code{lng}}{Present day longitude of the occurrence.}
#' 	\item{\code{paleolat}}{Reconstructed paleolatitude of the occurrence.}
#' 	\item{\code{paleolng}}{Reconstructed paleolongitude of the occurrence.}
#' 	\item{\code{period}}{Period of origin.}
#' 	\item{\code{epoch}}{Epoch of origin.}
#' 	\item{\code{subepoch}}{Subepoch of origin.}
#' 	\item{\code{stage}}{Geologic stage of the embedding rocks.}
#' 	\item{\code{early_interval}}{Early interval name registered in the PaleoDB dynamic time scale.}
#' 	\item{\code{late_interval}}{Late interval name registered in the PaleoDB dynamic time scale.}
#' 	\item{\code{max_ma}}{Maximum estimated age based on the PaleoDB dynamic time scale.}
#' 	\item{\code{min_ma}}{Minimum estimated age based on the PaleoDB dynamic time scale.}
#' 	\item{\code{slc}}{Bin number in the stage-level timescale \code{\link{stages}}.}
#' 	\item{\code{Bin}}{Bin number in the PaleoDB 10 million year resolution timescale \code{\link{bins}}.}
#' 	\item{\code{env}}{Environment of the occurrence: reefal \code{(r)}, non-reefal \code{(nr)} or unknown (\code{uk}), based on \code{\link{keys}}.}
#' 	\item{\code{lith}}{Substrate of the occurrence: carbonate \code{(c)}, siliciclastic \code{(s)} or unknown (\code{uk}), based on \code{\link{keys}}.}
#' 	\item{\code{latgroup}}{Latitude of the occurrence: tropical \code{(t)} or non-tropical \code{(nt)}.}
#' 	\item{\code{bath}}{Inferred depth of the occurrence: deep \code{(deep)}, shallow \code{(shal)} or unknown (\code{uk}), based on \code{\link{keys}}.}
#' 	\item{\code{gensp}}{The binomen of the occurrence.}
#' 	\item{\code{ecology}}{Symbiotic status of the occurring coral: zooxanthellate \code{(z)} or azooxanthellate \code{(az}, including apozooxanthellates).}
#' 	\item{\code{ecologyMostZ}}{Symbiotic status of the occurring coral, incorporating the uncertainty of inferred symbiotic status. This variable includes assignment with the maximum likely number of zooxanthellate genera.}
#' 	\item{\code{ecologyMostAZ}}{Symbiotic status of the occurring coral, incorporating the uncertainty of inferred symbiotic status. This variable includes assignment with the maximum likely number of azooxanthellate genera.}
#' 	\item{\code{growth}}{Growth type of the coral: \code{colonial} or \code{solitary}.}
#' 	\item{\code{integration}}{Integration of corallites from the scale of 0 to 4. \code{solitary} corals are marked with 0s.}
#' }
#' @source \url{http://www.paleobiodb.org/}
#' 
"corals"

#' Keys to process stratigraphic, environmental and lithological information from the Paleobiology Database
#' 
#' Lists of entries treated as indicators of similar characteristics
#' 
#' Entries in the stratigraphic, lithological and environment fields of current Paleobiology Database downloads are too numerous to form the basis of analyses without transformations. 
#' This variable includes potential groupings of entries that represent similar characteristics. These objects can be used by the \code{\link{map}} function to create new variables of stratigraphic, environmental and lithological information.
#' 
#'  
#' @format A \code{list} of 7 \code{list}s:
#' 	\describe{
#' 		\item{\code{binInt}}{A \code{list} of \code{vector}s. Entries in the \code{early_interval} and \code{late_interval} variables of PaleoDB downloads indicate the collections' positions in the dynamic time scale. These entries were linked to 10 million year-resolution time scale stored in \code{\link{bins}}. These links were compiled using a download from the FossilWorks website (\url{http://www.fossilworks.com/}), on 08 June, 2018. You can check the lookup table \code{\link{stratkeys}} here. This is version 0.9.1}
#' 		\item{\code{slcInt}}{A \code{list} of \code{vector}s. Entries in the \code{early_interval} and \code{late_interval} variables of PaleoDB downloads indicate the collections' positions in the dynamic time scale. These entries were linked to stage-resolution time scale stored in \code{\link{stages}}. See \code{binInt} for version information.} These entries are reliable only in the Post-Ordovician!
#' 		\item{\code{reefs}}{A \code{list} of \code{vector}s. Entries in the \code{environment} field of the PaleoDB download indicate information regarding the likely reefal origin of carbonatic rocks. See the vignette ('§PhaneroCurve') on the exact use of these data. v0.9.}
#' 		\item{\code{lith}}{A \code{list} of \code{vector}s. Entries in the \code{lithology1} field of the PaleoDB download indicate information regarding the substrate of the embedding rocks. This key maps the entries to \code{siliciclastic}, \code{"carbonate"} or \code{"unknown"} substrates. v0.9.}
#' 		\item{\code{lat}}{A \code{list} of \code{vector}s. Entries in the \code{paleolat} field of the PaleoDB download indicate information regarding paleolatitude of the occurrences. This key maps the entries to \code{"tropical"} or \code{"non-tropical"} latitudes. v0.9.}
#' 		\item{\code{grain}}{A \code{list} of \code{vector}s. Entries in the \code{lithology1} field of the PaleoDB download indicate information regarding the grain sizes of the depositional environment. This key maps the entries to \code{"coarse"}, \code{"fine"} or \code{"unknonwn"} grain sizes. v0.9.}
#' 
#' 	}
#' 
#' @source Stratigraphic assignments are based on the download of collection data from Fossilworks (\url{http://www.fossilworks.com/}) and the dynamic time scale of the Paleobiology Database, written by J. Alroy. The assignment of numeric values were done by A. Kocsis. Environmental variables were grouped by W. Kiessling.
"keys"

#' 95 bin Phanerozoic time scale based on the stratigraphic stages of Gradstein et al. (2012).
#' 
#' Stage-level timescale used in some analyses.
#' 
#' This is an example time scale object that can be used in the Phanerozoic scale analyses. Example occurrence datasets related to the package use the variable \code{slc} when referring to this timescale.
#' 
#' @format A \code{data.frame} with 95 observations and 9 variables:
#' 	\describe{
#' 		\item{\code{per}}{Abbreviation of geologic periods.}
#' 		\item{\code{period}}{Geologic periods.}
#' 		\item{\code{series}}{Geologic series.}
#' 		\item{\code{stage}}{Names of geologic stages.}
#' 		\item{\code{short}}{Abbreviations of geologic stages.}
#' 		\item{\code{bottom}}{Numeric ages of the bottoms boundaries (earliest ages) of the bins.}
#' 		\item{\code{mid}}{Numeric age midpoints of the bins, the averages of \code{bottom} and \code{top}.}
#' 		\item{\code{top}}{Numeric age of the tops (latest ages) of the bins.}
#' 		\item{\code{num}}{Integer number identifiers of the bins.}
#' 	}
#' 
#' @source Based on Gradstein et al. (2012), compiled by Wolfgang Kiessling.
"stages"


#' The 10 million year resolution timescale of the Paleobiology Database
#' 
#' Roughly 10 million year timescale used in some analyses. 
#' 
#' This is an example time scale object that can be used in the Phanerozoic scale analyses. This time scale comprises 49 bins, roughly 10 million years of durations that result from the combination of certain standard stages.
#' 
#' @format A \code{data.frame} with 95 observations and 9 variables:
#' 	\describe{
#' 		\item{X10}{The name of the bin: Period and number.}
#' 		\item{ocean}{The primary state of the oceans from the point of carbonate precipitation. \code{ar} indicates aragonitic, \code{cc} indicates calcitic conditions. Based on §}
#' 		\item{ocean2}{§}
#' 		\item{climate}{Primary climatic characteristic: \code{w} denotes warm, \code{c} denotes cold.}
#' 		\item{\code{bottom}}{Numeric ages of the bottom boundaries (earliest ages) of the bins.}
#'  		\item{\code{mid}}{Numeric ages midpoints of the bins, the averages of \code{bottom} and \code{top}.}
#'  		\item{\code{top}}{Numeric ages of the tops (latest ages) of the bins.}
#'  		\item{\code{BIN}}{Integer number identifiers of the bins. §correct to num!}
#' }
#' 
#' 
#' @source § executive committee of old PaleoDB, John Alroy. Additional variables were added by Wolfgang Kiessling.
"bins"

#' The FossilWorks-based lookup table for the stratigraphic assignments of collections in the Paleobiology Database
#' 
#' Table including the user-chosen interval data and the stratigraphic units of the dynamic timescale.
#' 
#' Since the separation of the FossilWorks (\url{http://www.fossilworks.com/}) portal from the Paleobiology Database (\url{http://www.paleobiodb.org/}) the access to the stratigraphic information in the database have been problematic. This table includes groupings of 
#' \code{early_interval}/\code{max_interval} entries of the dynamic timescale that users can choose during collection entry. The table assigns these intervals to some corresponding stratigraphic units from different time scales.
#' These entries were distilled from those collections that only have a \code{max_interval} value. As there is a mismatch between the data Paleobiology Database and FossilWorks this list is not comprehensive and a couple entries are probably missing. For this reason, this dataset is expected to be updated in the future. 
#' 
#' This particular version (v0.9.1) is based on a download of all collections in FossilWorks between the Ediacaran and the Holocene. The download took place on 22 June, 2018. The entries were transformed to \code{\link{keys}} to be used with the \code{\link{map}} function.
#' 
#' @format A \code{data.frame} with 758 observations of 7 variables:
#' 	\describe{
#' 		\item{\code{interval}}{The names of the registered intervals in the \code{early_interval}/\code{max_interval} and \code{late_interval}/\code{min_interval} columns.}
#' 		\item{\code{period}}{The period containing the interval.}
#' 		\item{\code{epoch}}{The epoch containing the interval.}
#' 		\item{\code{X10_my_bin}}{The 10 million year time scale interval containing the interval.}
#' 		\item{\code{bin}}{Numeric identifier of the 10 million year interval in the \code{\link{bins}} object.}
#' 		\item{\code{stage}}{The stage containing the interval.}	
#' 		\item{\code{slc}}{Numeric identifier of the interval in the stage-level time scale provided as \code{\link{stages}} object.}
#' 	}
#' 	
#' @source \url{http://www.fossilworks.com/}
"stratkeys"