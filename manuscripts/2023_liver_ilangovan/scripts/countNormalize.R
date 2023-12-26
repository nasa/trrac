# library imports
library(DESeq2)
library(dplyr)
library(optparse)

# the arguments provided at runtime for 
# feature subset (F), metadata directory (M), raw counts directory (C)
option_list <- list(
    make_option(c("-F", "--feature_list"),
    type = "character",
    default = NULL,
    help = "gene set filtering file name",
    metavar = "character"),
    make_option(c("-M", "--metadata_dir"),
    type = "character",
    default = "single_study_metadata",
    help = "metadata dir name",
    metavar = "character"),
    make_option(c("-C", "--counts_dir"),
    type = "character",
    default = "raw_counts",
    help = "metadata dir name",
    metavar = "character")
); 

# parses the command line options provided at script runtime
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# logging the filepaths used for feature subset, metadata, and counts
if(!is.null(opt$feature_list)) {
    print(paste(opt$feature_list, opt$metadata_dir, opt$counts_dir));
};

# subsets the counts matrix based on the provided feature list
parse_features <- function(counts) {
    if (!is.null(opt$feature_list)) {
        features <- read.table(opt$feature_list, header = TRUE, sep = ",")
        print(paste(features[1,2]))
        counts <- counts[row.names(counts) %in% features[, 2], ]
    };

    return(counts)
}

# loads the counts and metadata files
file_imports <- function(counts_fname, metadata_fname) {
    counts_path <- paste(getwd(),
                        "data",
                        "raw_counts",
                         paste0(counts_fname, ".csv"),
                         sep = "/")

    metadata_path <- paste(getwd(),
                            "data",
                            "single_study_metadata",
                            paste0(metadata_fname, ".txt"),
                            sep = "/")

    counts <- read.delim(counts_path,
                            sep = ",",
                            row.names = 1,
                            header = TRUE)

    coldata <- read.delim(metadata_path,
                            row.names = 1,
                            header = TRUE,
                            sep = "\t")

    return(list(counts = counts, coldata = coldata))
}

### data set creation function for each GLDS

# filtering and counts loading
glds47_dds <- function(counts_fname, metadata_fname) {
    data <- file_imports(counts_fname, metadata_fname)

    counts <- data$counts
    coldata <- data$coldata

    counts <- parse_features(counts)

    # clean the factor value names by replacing " " -> "_" and removing "#"
    coldata$spaceflight <- coldata$Factor.Value.Spaceflight. %>%
                            { gsub(" ", "_", .) } %>% # nolint
                            { gsub("#", "", .) } # nolint

    # remove Basal Control samples from the DGE analysis
    coldata <- coldata[!coldata$spaceflight %in% "Basal_Control", ]
    coldata$Sample.Name <- coldata$Sample.Name %>% # nolint
                            { gsub("-", ".", .) } %>% # nolint
                            factor()

    counts <- counts[, colnames(counts) %in% coldata$Sample.Name]

    coldata <- coldata[match(colnames(counts), coldata$Sample.Name), ]

    # apply factor conversion after filtering to omit filtered out factor levels
    coldata$spaceflight <- factor(coldata$spaceflight)

    # dimensions of data after filtering out basal control
    cat("Dimensions of Counts:", dim(counts), "\n")
    cat("Dimensions of Metadata:", dim(coldata), "\n")

    # Create dataset for DESeq2
    # design factor is spaceflight
    factor_names <- c("spaceflight")
    design_formula <- paste("",
                            paste(factor_names,
                            collapse = " + "),
                            sep = " ~ ")

    cat("The design formula is:", design_formula)

    row.names(coldata) <- coldata$Sample.Name

    dds_nofilt <- DESeqDataSetFromMatrix(countData = round(counts),
                                        colData = coldata,
                                        design = as.formula(design_formula))

    return(list(data = dds_nofilt, metadata = coldata))
}

# filtering and counts loading
glds242_dds <- function(counts_fname, metadata_fname) {
    data <- file_imports(counts_fname, metadata_fname)

    counts <- data$counts
    coldata <- data$coldata

    # remove ERCC spikes
    ercc_spikes <- row.names(counts) %>%
                    { grep("^ERCC", .) } #nolint
    counts <- counts[-ercc_spikes, ]
    
    counts <- parse_features(counts)

    # clean the factor value names by replacing " " -> "_" and removing "#"
    coldata$spaceflight <- coldata$Factor.Value.Spaceflight. %>%
                            { gsub(" ", "_", .) } %>% #nolint
                            { gsub("#", "", .) } #nolint

    # remove Basal Control and Vivarium samples from DGE analysis
    coldata <- coldata[!coldata$spaceflight %in%
                c("Basal_Control", "Vivarium_Control"), ]

    coldata$Sample.Name <- coldata$Sample.Name %>%
                            { gsub("-", ".", .) } %>% #nolint
                            factor()
    counts <- counts[,colnames(counts) %in% coldata$Sample.Name]
    coldata <- coldata[match(coldata$Sample.Name, colnames(counts)), ]

    # apply factor conversion after filtering to omit filtered out factor levels
    coldata$spaceflight <- factor(coldata$spaceflight) %>%
                            relevel(ref = "Ground_Control")

    # dimensions of data after filtering out basal control
    cat("Dimensions of Counts:", dim(counts), "\n")
    cat("Dimensions of Metadata:", dim(coldata), "\n")

    # Create dataset for DESeq2
    # design factor is spaceflight
    factor_names <- c("spaceflight")
    design_formula <- paste("",
                        paste(factor_names, collapse = " + "),
                        sep = " ~ ")

    cat("The design formula is:", design_formula)

    row.names(coldata) <- coldata$Sample.Name

    dds_nofilt <- DESeqDataSetFromMatrix(countData = round(counts),
                                        colData = coldata,
                                        design = as.formula(design_formula))

    return(list(data = dds_nofilt, metadata = coldata))
}

# filtering and counts loading
glds245_dds <- function(counts_fname, metadata_fname) {
    data <- file_imports(counts_fname, metadata_fname)

    counts <- data$counts
    coldata <- data$coldata

    extract_path <- paste(getwd(),
                        "data",
                        "single_study_metadata",
                        "specialty",
                         paste0("245_RR6_extraction dates", ".txt"),
                         sep = "/")

    extractdata <- read.delim(extract_path,
                            row.names = 1,
                            header = TRUE,
                            sep = "\t")

    # remove ERCC spikes
    ercc_spikes <- row.names(counts) %>%
                    { grep("^ERCC", .) } #nolint
    counts <- counts[-ercc_spikes, ]

    counts <- parse_features(counts)

    # clean the factor value names by replacing " " -> "_" and removing "#"
    coldata$spaceflight <- coldata$Factor.Value.Spaceflight. %>%
                            { gsub(" ", "_", .) } %>% #nolint
                            { gsub("#", "", .) } #nolint
    coldata$duration <- coldata$Factor.Value.Duration. %>%
                        { gsub("~", "", .) } #nolint
    coldata$euthanasia <- coldata$Factor.Value.Euthanasia. %>%
                            { gsub(" ", "_", .) } #nolint
    coldata$dissection <- coldata$Factor.Value.Dissection.Condition. %>%
                            { gsub(" ", "_", .) } #nolint

    # remove Basal Control samples from the DGE analysis
    coldata <- coldata[!coldata$spaceflight %in% "Basal_Control",]
    outliers <- c("LAR Flight 5", "ISS-T Flight 5", "ISS-T Flight 9")
    coldata <- coldata[!(rownames(coldata) %in% outliers), ]

    # adjusting for batch effect
    extractdata <- extractdata[rownames(extractdata) %in% coldata$Sample.Name,]

    coldata$Sample.Name <- coldata$Sample.Name %>%
                            { gsub("-", ".", .) } %>% #nolint
                            factor()
    counts <- counts[,colnames(counts) %in% coldata$Sample.Name]
    coldata <- coldata[match(coldata$Sample.Name, colnames(counts)), ]
    coldata$month_year <- extractdata$Month.year

    # apply factor conversion after filtering to omit filtered out factor levels
    coldata$spaceflight <- factor(coldata$spaceflight)
    coldata$duration <- factor(coldata$duration)
    coldata$euthanasia <- factor(coldata$euthanasia)
    coldata$dissection <- factor(coldata$dissection)
    coldata$month_year <- factor(coldata$month_year)

    # dimensions of data after filtering out basal control
    cat("Dimensions of Counts:", dim(counts), "\n")
    cat("Dimensions of Metadata:", dim(coldata), "\n")

    # Create dataset for DESeq2
    # cannot include all factors without design matrix losing full rank
    factor_names <- c("spaceflight", "duration", "euthanasia", "month_year")
    design_formula <- paste("",
                            paste(factor_names, collapse = " + "),
                            sep = " ~ ")

    cat("The design formula is:", design_formula)

    row.names(coldata) <- coldata$Sample.Name

    dds_nofilt <- DESeqDataSetFromMatrix(countData = round(counts),
                                        colData = coldata,
                                        design = as.formula(design_formula))

    return(list(data = dds_nofilt, metadata = coldata))
}

glds168_rr1_dds <- function(counts_fname, metadata_fname) {
    counts_path <- paste(getwd(),
                        "data",
                        "raw_counts",
                         paste0(counts_fname, ".csv"),
                         sep = "/")

    metadata_path <- paste(getwd(),
                            "data",
                            "single_study_metadata",
                            "specialty",
                            paste0(metadata_fname, ".txt"),
                            sep = "/")

    counts <- read.delim(counts_path,
                            sep = ",",
                            row.names = 1,
                            header = TRUE)

    coldata <- read.delim(metadata_path,
                            row.names = 2,
                            header = TRUE,
                            sep = "\t")

    # remove ERCC spikes
    ercc_spikes <- row.names(counts) %>%
                    { grep("^ERCC", .) } #nolint
    counts <- counts[-ercc_spikes, ]

    counts <- parse_features(counts)

    # clean the factor value names by replacing " " -> "_" and removing "#"
    coldata$spaceflight <- coldata$Factor.Value.Spaceflight. %>%
                            { gsub(" ", "_", .) } %>% #nolint
                            { gsub("#", "", .) } #nolint

    coldata$mission <- coldata$Factor.Value.Space.Mission. %>%
                        { gsub(" ", "_", .) } #nolint

    coldata$spikein <- coldata$Factor.Value.Spike.in.Quality.Control. %>%
                        { gsub(" ", "_", .) } #nolint

    # remove the Basal Control samples from the DGE analysis
    coldata <- coldata[!coldata$spaceflight %in%
                        c("Basal_Control", "Vivarium_Control"), ]

    coldata$Sample.Name <- row.names(coldata) %>%
                            { gsub("-", ".", .) } %>% #nolint
                            factor()

    counts<- counts[,colnames(counts) %in%
                    coldata$Sample.Name]

    coldata <- coldata[match(coldata$Sample.Name,
                            colnames(counts)), ]

    # apply factor conversion after filtering to omit filtered out factor levels
    coldata$spaceflight <- factor(coldata$spaceflight)
    coldata$mission <- factor(coldata$mission)
    coldata$spikein <- factor(coldata$spikein)

    missions <- c("SpaceX-4_(RR1)", "SpaceX-8_(RR3)")

    rr1_coldata <- coldata %>%
                    dplyr::filter(mission == "SpaceX-4_(RR1)")

    rr1_counts <- counts[,colnames(counts) %in%
                            rr1_coldata$Sample.Name]

    rr1_coldata <- rr1_coldata[match(colnames(rr1_counts),
                                rr1_coldata$Sample.Name),]

    counts_ercc_spikes <- colnames(rr1_counts) %>%
                    { grep("noERCC", .) } #nolint
    cols_ercc_spikes <- rr1_coldata$Sample.Name %>%
                    { grep("noERCC", .) } #nolint
    rr1_counts <- rr1_counts[,-counts_ercc_spikes]
    rr1_coldata <- rr1_coldata[-cols_ercc_spikes,]

    # dimensions of data after filtering out basal control
    cat("Dimensions of RR1 Counts:", dim(rr1_counts), "\n")
    cat("Dimensions of RR1 Metadata:", dim(rr1_coldata), "\n")

    # Create dataset for DESeq2
    # cannot include all factors without design matrix losing full rank
    all_factor_names = c("spaceflight", "mission", "spikein")
    factor_names = c("spaceflight")
    design_formula <- paste("", 
                        paste(factor_names,
                            collapse = " + "),
                            sep = " ~ ")

    cat("The design formula is:", design_formula)

    row.names(rr1_coldata) <- rr1_coldata$Sample.Name
    
    dds_nofilt <- DESeqDataSetFromMatrix(countData = round(rr1_counts),
                                            colData = rr1_coldata,
                                            design = as.formula(design_formula))

    return(list(data = dds_nofilt, metadata = rr1_coldata))
}

# filtering and counts loading
glds168_rr3_dds <- function(counts_fname, metadata_fname) {
    counts_path <- paste(getwd(),
                        "data",
                        "raw_counts",
                         paste0(counts_fname, ".csv"),
                         sep = "/")

    metadata_path <- paste(getwd(),
                            "data",
                            "single_study_metadata",
                            "specialty",
                            paste0(metadata_fname, ".txt"),
                            sep = "/")

    counts <- read.delim(counts_path,
                            sep = ",",
                            row.names = 1,
                            header = TRUE)

    coldata <- read.delim(metadata_path,
                            row.names = 2,
                            header = TRUE,
                            sep = "\t")

    # remove ERCC spikes
    ercc_spikes <- row.names(counts) %>%
                    { grep("^ERCC", .) } #nolint
    counts <- counts[-ercc_spikes, ]

    counts <- parse_features(counts)

    # clean the factor value names by replacing " " -> "_" and removing "#"
    coldata$spaceflight <- coldata$Factor.Value.Spaceflight. %>%
                            { gsub(" ", "_", .) } %>% #nolint
                            { gsub("#", "", .) } #nolint

    coldata$mission <- coldata$Factor.Value.Space.Mission. %>%
                        { gsub(" ", "_", .) } #nolint

    coldata$spikein <- coldata$Factor.Value.Spike.in.Quality.Control. %>%
                        { gsub(" ", "_", .) } #nolint

    # remove the Basal Control samples from the DGE analysis
    coldata <- coldata[!coldata$spaceflight %in%
                        c("Basal_Control", "Vivarium_Control"), ]

    coldata$Sample.Name <- row.names(coldata) %>%
                            { gsub("-", ".", .) } %>% #nolint
                            factor()

    counts<- counts[,colnames(counts) %in%
                    coldata$Sample.Name]

    coldata <- coldata[match(coldata$Sample.Name,
                            colnames(counts)), ]

    # apply factor conversion after filtering to omit filtered out factor levels
    coldata$spaceflight <- factor(coldata$spaceflight)
    coldata$mission <- factor(coldata$mission)
    coldata$spikein <- factor(coldata$spikein)

    missions <- c("SpaceX-4_(RR1)", "SpaceX-8_(RR3)")

    rr3_coldata <- coldata %>%
                    dplyr::filter(mission == "SpaceX-8_(RR3)")

    rr3_counts <- counts[,colnames(counts) %in%
                            rr3_coldata$Sample.Name]

    rr3_coldata <- rr3_coldata[match(colnames(rr3_counts),
                                rr3_coldata$Sample.Name),]

    # dimensions of data after filtering out basal control
    cat("Dimensions of RR3 Counts:", dim(rr3_counts), "\n")
    cat("Dimensions of RR3 Metadata:", dim(rr3_coldata), "\n")

    # Create dataset for DESeq2
    # cannot include all factors without design matrix losing full rank
    all_factor_names = c("spaceflight", "mission", "spikein")
    factor_names = c("spaceflight")
    design_formula <- paste("",
                        paste(factor_names,
                            collapse = " + "),
                            sep = " ~ ")

    cat("The design formula is:", design_formula)

    row.names(rr3_coldata) <- rr3_coldata$Sample.Name
    
    dds_nofilt <- DESeqDataSetFromMatrix(countData = round(rr3_counts),
                                            colData = rr3_coldata,
                                            design = as.formula(design_formula))

    return(list(data = dds_nofilt, metadata = rr3_coldata))
}

# filtering and counts loading
glds379_dds <- function(counts_fname, metadata_fname) {
    data <- file_imports(counts_fname, metadata_fname)

    counts <- data$counts
    coldata <- data$coldata

    # remove ERCC spikes
    ercc_spikes <- row.names(counts) %>%
                    { grep("^ERCC", .) } #nolint
    counts <- counts[-ercc_spikes, ]

    counts <- parse_features(counts)

    # clean the factor value names by replacing " " -> "_" and removing "#"
    coldata$spaceflight <- coldata$Factor.Value.Spaceflight. %>%
                            { gsub(" ", "_", .) } %>% #nolint
                            { gsub("#", "", .) } #nolint

    # remove Basal Control samples from the DGE analysis
    coldata <- coldata[!coldata$spaceflight %in% c("Vivarium_Control","Basal_Control"),]
    outliers <- c('FL-ISS-10','FL-ISS-16','FL-ISS-17','FL-LAR-02',
              'FL-LAR-08','FL-LAR-11','FL-LAR-15','FL-LAR-18',
              'HGC-LAR-03','HGC-LAR-05','HGC-LAR-18')
    coldata <- coldata[!(rownames(coldata) %in% outliers), ]

    coldata$Sample.Name <- coldata$Sample.Name %>%
                            { gsub("-", ".", .) } 

    counts <- counts[,colnames(counts) %in% coldata$Sample.Name]
    coldata <- coldata[match(coldata$Sample.Name, colnames(counts)), ]

    # apply factor conversion after filtering to omit filtered out factor levels
    coldata$spaceflight <- factor(coldata$spaceflight)

    # dimensions of data after filtering out basal control
    cat("Dimensions of Counts:", dim(counts), "\n")
    cat("Dimensions of Metadata:", dim(coldata), "\n")

    # Create dataset for DESeq2
    # cannot include all factors without design matrix losing full rank
    factor_names <- c("spaceflight")
    design_formula <- paste("",
                            paste(factor_names, collapse = " + "),
                            sep = " ~ ")

    cat("The design formula is:", design_formula)

    row.names(coldata) <- coldata$Sample.Name

    dds_nofilt <- DESeqDataSetFromMatrix(countData = round(counts),
                                        colData = coldata,
                                        design = as.formula(design_formula))

    return(list(data = dds_nofilt, metadata = coldata))
}

# generate DESeq objects from counts and metadata
glds47_output <- glds47_dds("47", "47")
glds242_output <- glds242_dds("242", "242")
glds245_output <- glds245_dds("245", "245")
glds168_rr1_output <- glds168_rr1_dds("168", "168")
glds168_rr3_output <- glds168_rr3_dds("168", "168")
glds379_output <- glds379_dds("379", "379")

# perform factor size estimation (median of ratios normalization)
glds47_dds <- estimateSizeFactors(glds47_output$data)
glds242_dds <- estimateSizeFactors(glds242_output$data)
glds245_dds <- estimateSizeFactors(glds245_output$data)
glds168_rr1_dds <- estimateSizeFactors(glds168_rr1_output$data)
glds168_rr3_dds <- estimateSizeFactors(glds168_rr3_output$data)
glds379_dds <- estimateSizeFactors(glds379_output$data)

# extract the normalized counts files
norm_47 <- counts(glds47_dds, normalized = TRUE)
norm_242 <- counts(glds242_dds, normalized = TRUE)
norm_245 <- counts(glds245_dds, normalized = TRUE)
norm_168_rr1 <- counts(glds168_rr1_dds, normalized = TRUE)
norm_168_rr3 <- counts(glds168_rr3_dds, normalized = TRUE)
norm_379 <- counts(glds379_dds, normalized = TRUE)

colnames(norm_47) <- colnames(norm_47) %>%
                        { gsub("\\.", "-", .) } #nolint
colnames(norm_242) <- colnames(norm_242) %>%
                        { gsub("\\.", "-", .) } #nolint
colnames(norm_245) <- colnames(norm_245) %>%
                        { gsub("\\.", "-", .) } #nolint
colnames(norm_168_rr1) <- colnames(norm_168_rr1) %>%
                        { gsub("\\.", "-", .) } #nolint                        
colnames(norm_168_rr3) <- colnames(norm_168_rr3) %>%
                        { gsub("\\.", "-", .) } #nolint     
colnames(norm_379) <- colnames(norm_379) %>%
                        { gsub("\\.", "-", .) } #nolint                        

# write out files into the normalized counts path
norm_path <- paste(getwd(),
                        "data",
                        "norm_counts",
                         sep = "/")

write.table(round(norm_47),
            file = paste0(norm_path,
                    "/",
                    "47.csv"),
            sep = ",",
            row.names = TRUE,
            col.names = NA)

write.table(round(norm_242),
            file = paste0(norm_path,
                    "/",
                    "242.csv"),
            sep = ",",
            row.names = TRUE,
            col.names = NA)

write.table(round(norm_245),
            file = paste0(norm_path,
                    "/",
                    "245.csv"),
            sep = ",",
            row.names = TRUE,
            col.names = NA)

write.table(round(norm_168_rr1),
            file = paste0(norm_path,
                    "/",
                    "168_rr1.csv"),
            sep = ",",
            row.names = TRUE,
            col.names = NA)

write.table(round(norm_168_rr3),
            file = paste0(norm_path,
                    "/",
                    "168_rr3.csv"),
            sep = ",",
            row.names = TRUE,
            col.names = NA)

write.table(round(norm_379),
            file = paste0(norm_path,
                    "/",
                    "379.csv"),
            sep = ",",
            row.names = TRUE,
            col.names = NA)
