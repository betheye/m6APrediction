#' Encode DNA Sequences into Numerical Features
#'
#' @description
#' This function converts DNA sequences into a data frame with each nucleotide
#' position as a separate factor column. It is used internally for feature
#' engineering in the m6A prediction pipeline.
#'
#' @param dna_strings A character vector containing DNA sequences of equal length
#'
#' @return A data frame where each column represents a nucleotide position,
#'   with factor levels A, T, C, G
#'
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Multiple Sample Prediction for m6A Sites
#'
#' @description
#' This function predicts m6A modification sites for multiple samples using
#' a pre-trained random forest model. It takes a data frame of features and
#' returns predictions with probability scores.
#'
#' @param ml_fit A trained random forest model object (e.g., from randomForest package)
#' @param feature_df A data frame containing the following required columns:
#'   \describe{
#'     \item{gc_content}{Numeric, GC content proportion (0-1)}
#'     \item{RNA_type}{Character, one of "mRNA", "lincRNA", "lncRNA", "pseudogene"}
#'     \item{RNA_region}{Character, one of "CDS", "intron", "3'UTR", "5'UTR"}
#'     \item{exon_length}{Numeric, length of exon in nucleotides}
#'     \item{distance_to_junction}{Numeric, distance to nearest splice junction}
#'     \item{evolutionary_conservation}{Numeric, conservation score (0-1)}
#'     \item{DNA_5mer}{Character, 5-nucleotide DNA sequence (A/T/C/G)}
#'   }
#' @param positive_threshold A numeric value between 0 and 1 for classification
#'   threshold. Sites with probability above this value are classified as
#'   "Positive" (default: 0.5)
#'
#' @return A data frame with two additional columns: predicted_m6A_prob
#'   (the predicted probability) and predicted_m6A_status (either "Positive"
#'   or "Negative" based on the threshold)
#' @importFrom stats predict
#' @import randomForest
#' @export
#'
#' @examples
#' # Load the pre-trained random forest model
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'
#' # Load example input data
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                     package = "m6APrediction"))
#'
#' # Make predictions with default threshold (0.5)
#' predictions <- prediction_multiple(ml_fit, example_df)
#'
#' # Use a more stringent threshold
#' predictions_stringent <- prediction_multiple(ml_fit, example_df, positive_threshold = 0.7)
#'
#' # View the first few predictions
#' head(predictions)
#'
#' @seealso \code{\link{prediction_single}} for single sample prediction
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length",
                  "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
                %in% colnames(feature_df)))

  if (!inherits(ml_fit, "randomForest")) {
    stop("ml_fit must be a randomForest model object")
  }

  if (!is.numeric(positive_threshold) || positive_threshold < 0 || positive_threshold > 1) {
    stop("positive_threshold must be a numeric value between 0 and 1")
  }

  if (nrow(feature_df) == 0) {
    stop("feature_df must contain at least one row")
  }

  encoded_seq <- dna_encoding(feature_df$DNA_5mer)
  pred_df <- cbind(feature_df, encoded_seq)

  pred_df$RNA_type <- factor(pred_df$RNA_type,
                             levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  pred_df$RNA_region <- factor(pred_df$RNA_region,
                               levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  pred_probs <- predict(ml_fit, newdata = pred_df, type = "prob")
  feature_df$predicted_m6A_prob <- round(pred_probs[, "Positive"], 4)
  feature_df$predicted_m6A_status <- ifelse(
    feature_df$predicted_m6A_prob > positive_threshold,
    "Positive",
    "Negative"
  )

  return(feature_df)
}


#' Single Sample Prediction for m6A Sites
#'
#' @description
#' This function predicts m6A modification probability and status for a single
#' sample based on individual feature values. It validates the DNA sequence
#' input and calls the prediction_multiple function internally.
#'
#' @param ml_fit A trained random forest model object
#' @param gc_content A numeric value representing the GC content
#' @param RNA_type A character string specifying the RNA type (one of: "mRNA",
#'   "lincRNA", "lncRNA", "pseudogene")
#' @param RNA_region A character string specifying the RNA region (one of:
#'   "CDS", "intron", "3'UTR", "5'UTR")
#' @param exon_length A numeric value for the exon length
#' @param distance_to_junction A numeric value for the distance to the
#'   nearest junction
#' @param evolutionary_conservation A numeric value representing the
#'   conservation score
#' @param DNA_5mer A character string of exactly 5 nucleotides (A, T, C, G)
#'   representing the 5-mer DNA sequence
#' @param positive_threshold A numeric value between 0 and 1 for classification
#'   threshold (default: 0.5)
#'
#' @return A named vector with two elements: predicted_m6A_prob (the predicted
#'   probability) and predicted_m6A_status (either "Positive" or "Negative")
#'
#' @export
#'
#' @examples
#' # Load the pre-trained random forest model
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'
#' # Make a prediction for a single sample
#' result <- prediction_single(
#'   ml_fit = ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "ATCGA",
#'   positive_threshold = 0.5
#' )
#'
#' # View the prediction result
#' print(result)
#'
#' @seealso \code{\link{prediction_multiple}} for multiple sample prediction
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5){

  # Validate DNA sequence
  if(nchar(DNA_5mer) != 5){
    stop("DNA_5mer must be exactly 5 characters long")
  }
  if(!grepl("^[ATCGatcg]+$", DNA_5mer)){
    stop("DNA_5mer must contain only A, T, C, G characters")
  }

  DNA_5mer <- toupper(DNA_5mer)

  single_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  result_df <- prediction_multiple(ml_fit, single_df, positive_threshold)

  returned_vector <- c(
    predicted_m6A_prob = result_df$predicted_m6A_prob[1],
    predicted_m6A_status = result_df$predicted_m6A_status[1]
  )

  return(returned_vector)
}
