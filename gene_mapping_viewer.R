#!/usr/bin/env Rscript

# Gene Mapping Visualization Shiny App
# Visualizes gene mappings between ICSASG_v2 and Ssal_v3.1 assemblies

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Function to check if a seqid is a chromosome (not a scaffold)
is_chromosome <- function(seqid) {
  # Match numeric chromosomes (e.g., "1", "2", "3") or ssa-prefixed (e.g., "ssa01", "ssa1")
  # Excludes scaffolds, contigs, and unplaced sequences
  str_detect(seqid, "^(ssa)?\\d+$")
}

# Function to parse GFF3 and extract gene positions
parse_gff_genes <- function(gff_file) {
  message("Reading GFF file: ", gff_file)

  # Read GFF3 file, skip comment lines
  gff <- read_tsv(gff_file,
                  col_names = c("seqid", "source", "type", "start", "end",
                                "score", "strand", "phase", "attributes"),
                  comment = "#",
                  col_types = cols(.default = "c"))

  # Filter for gene features only
  genes <- gff %>%
    filter(type == "gene") %>%
    mutate(
      start = as.numeric(start),
      end = as.numeric(end),
      # Extract gene ID from attributes
      gene_id = str_extract(attributes, "ID=gene:([^;]+)") %>%
        str_replace("ID=gene:", "")
    ) %>%
    select(gene_id, seqid, start, end, strand)

  return(genes)
}

# Function to parse comparison TSV
parse_comparison <- function(tsv_file) {
  message("Reading comparison file: ", tsv_file)

  comp <- read_tsv(tsv_file, col_types = cols(.default = "c"))

  # Filter for gene-level comparisons only
  genes <- comp %>%
    filter(feature == "gene") %>%
    mutate(
      geneA = str_replace(geneA, "gene:", ""),
      geneB = str_replace(geneB, "gene:", ""),
      # Extract class from stats column
      class = str_extract(stats, "class=([^;]+)") %>%
        str_replace("class=", "")
    ) %>%
    select(geneA, geneB, class)

  return(genes)
}

# UI
ui <- fluidPage(
  titlePanel("Gene Mapping Viewer: ICSASG_v2 to Ssal_v3.1"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Data Files"),
      textInput("icsasg_gff", "ICSASG_v2 GFF:",
                value = "data/genomes/AtlanticSalmon/ICSASG_v2_Ens.gff3"),
      textInput("ssal_lifted_gff", "Lifted GFF:",
                value = "experiments/liftoff_full/ICSASG_v2_to_Ssal_v3.1_Ens.gff3"),
      textInput("ssal_native_gff", "Ssal_v3.1 Native GFF:",
                value = "data/genomes/AtlanticSalmon/Ssal_v3.1_Ens.gff3"),
      textInput("comparison_tsv", "Comparison TSV:",
                value = "experiments/comparison_runs/ens_lift_vs_native.tsv"),
      actionButton("load_data", "Load Data", class = "btn-primary"),
      hr(),
      textOutput("data_status"),
      hr(),
      h4("Filters"),
      checkboxInput("show_mapped", "Liftoff mapped (all)", value = TRUE),
      checkboxInput("show_id_match", "Liftoff mapped (ID match only)", value = FALSE),
      checkboxInput("show_stable_only", "Stable ID only (not in liftoff)", value = TRUE),
      helpText("Note: If both 'all' and 'ID match only' are checked, 'all' takes precedence."),
      hr(),
      h4("Selected Chromosome Pair"),
      textOutput("selected_pair"),
      textOutput("gene_count")
    ),

    mainPanel(
      width = 9,
      plotOutput("matrix_plot",
                 height = "600px",
                 click = "matrix_click"),
      hr(),
      plotOutput("dotplot", height = "600px")
    )
  )
)

# Server
server <- function(input, output, session) {

  # Reactive values to store data
  data <- reactiveValues(
    icsasg_genes = NULL,
    ssal_lifted_genes = NULL,
    ssal_native_genes = NULL,
    comparison = NULL,
    merged_data = NULL,
    matrix_data = NULL,
    selected_icsasg_chr = NULL,
    selected_ssal_chr = NULL
  )

  # Load data when button is clicked
  observeEvent(input$load_data, {
    withProgress(message = 'Loading data...', {

      incProgress(0.2, detail = "Reading ICSASG_v2 GFF")
      data$icsasg_genes <- parse_gff_genes(input$icsasg_gff)

      incProgress(0.4, detail = "Reading Ssal_v3.1 lifted GFF")
      data$ssal_lifted_genes <- parse_gff_genes(input$ssal_lifted_gff)

      incProgress(0.6, detail = "Reading Ssal_v3.1 native GFF")
      data$ssal_native_genes <- parse_gff_genes(input$ssal_native_gff)

      incProgress(0.8, detail = "Reading comparison data")
      data$comparison <- parse_comparison(input$comparison_tsv)

      incProgress(0.9, detail = "Merging data")

      # Merge comparison with position data (liftoff mappings)
      merged <- data$comparison %>%
        left_join(data$icsasg_genes, by = c("geneA" = "gene_id")) %>%
        rename(icsasg_chr = seqid, icsasg_start = start,
               icsasg_end = end, icsasg_strand = strand) %>%
        left_join(data$ssal_lifted_genes, by = c("geneA" = "gene_id")) %>%
        rename(ssal_lifted_chr = seqid, ssal_lifted_start = start,
               ssal_lifted_end = end, ssal_lifted_strand = strand) %>%
        left_join(data$ssal_native_genes, by = c("geneB" = "gene_id")) %>%
        rename(ssal_native_chr = seqid, ssal_native_start = start,
               ssal_native_end = end, ssal_native_strand = strand) %>%
        # Remove genes without position data
        filter(!is.na(icsasg_chr) & !is.na(ssal_lifted_chr)) %>%
        # Filter to only include chromosomes (not scaffolds)
        filter(is_chromosome(icsasg_chr) & is_chromosome(ssal_lifted_chr)) %>%
        # Add mapping type flags
        mutate(
          id_match = (geneA == geneB),
          mapping_type = "liftoff_mapped"
        )

      # Find stable ID only genes (same ID in both assemblies but not in comparison)
      icsasg_chr_genes <- data$icsasg_genes %>%
        filter(is_chromosome(seqid))

      ssal_native_chr_genes <- data$ssal_native_genes %>%
        filter(is_chromosome(seqid))

      # Genes with matching IDs in both assemblies
      stable_id_genes <- icsasg_chr_genes %>%
        inner_join(ssal_native_chr_genes, by = "gene_id", suffix = c("_icsasg", "_ssal")) %>%
        # Exclude genes already in comparison (by geneA)
        anti_join(merged, by = c("gene_id" = "geneA")) %>%
        mutate(
          geneA = gene_id,
          geneB = gene_id,
          class = "StableIDOnly",
          icsasg_chr = seqid_icsasg,
          icsasg_start = start_icsasg,
          icsasg_end = end_icsasg,
          icsasg_strand = strand_icsasg,
          ssal_lifted_chr = seqid_ssal,  # Use native position as "lifted"
          ssal_lifted_start = start_ssal,
          ssal_lifted_end = end_ssal,
          ssal_lifted_strand = strand_ssal,
          ssal_native_chr = seqid_ssal,
          ssal_native_start = start_ssal,
          ssal_native_end = end_ssal,
          ssal_native_strand = strand_ssal,
          id_match = TRUE,
          mapping_type = "stable_id_only"
        ) %>%
        select(geneA, geneB, class, icsasg_chr, icsasg_start, icsasg_end, icsasg_strand,
               ssal_lifted_chr, ssal_lifted_start, ssal_lifted_end, ssal_lifted_strand,
               ssal_native_chr, ssal_native_start, ssal_native_end, ssal_native_strand,
               id_match, mapping_type)

      # Combine both datasets
      merged_all <- bind_rows(merged, stable_id_genes) %>%
        # Calculate gene midpoints for plotting
        mutate(
          icsasg_pos = (icsasg_start + icsasg_end) / 2,
          ssal_pos = (ssal_lifted_start + ssal_lifted_end) / 2
        )

      data$merged_data <- merged_all

      incProgress(1.0, detail = "Done!")
    })
  })

  # Reactive filtered data based on checkboxes
  filtered_data <- reactive({
    req(data$merged_data)

    filtered <- data$merged_data

    # Apply filters based on checkboxes
    if (!input$show_mapped && !input$show_id_match && !input$show_stable_only) {
      # If all unchecked, return empty dataset
      return(filtered %>% filter(FALSE))
    }

    # Build filter conditions
    keep_rows <- rep(FALSE, nrow(filtered))

    if (input$show_mapped) {
      # Show all liftoff mapped genes (regardless of ID match)
      keep_rows <- keep_rows | (filtered$mapping_type == "liftoff_mapped")
    }

    if (input$show_id_match) {
      # Show only liftoff mapped genes where ID matches
      keep_rows <- keep_rows | (filtered$mapping_type == "liftoff_mapped" & filtered$id_match)
    }

    if (input$show_stable_only) {
      # Show stable ID only genes
      keep_rows <- keep_rows | (filtered$mapping_type == "stable_id_only")
    }

    filtered <- filtered[keep_rows, ]

    return(filtered)
  })

  # Data status text
  output$data_status <- renderText({
    if (is.null(data$merged_data)) {
      "No data loaded. Click 'Load Data' to begin."
    } else {
      total <- nrow(data$merged_data)
      filtered_count <- nrow(filtered_data())
      liftoff_count <- sum(data$merged_data$mapping_type == "liftoff_mapped")
      stable_count <- sum(data$merged_data$mapping_type == "stable_id_only")
      id_match_count <- sum(data$merged_data$mapping_type == "liftoff_mapped" & data$merged_data$id_match)

      paste0("Total: ", total, " (", liftoff_count, " liftoff, ",
             id_match_count, " ID match, ", stable_count, " stable ID only) | ",
             "Showing: ", filtered_count)
    }
  })

  # Matrix plot
  output$matrix_plot <- renderPlot({
    req(filtered_data())

    # Create matrix data for heatmap from filtered data
    matrix_data <- filtered_data() %>%
      count(icsasg_chr, ssal_lifted_chr) %>%
      complete(icsasg_chr, ssal_lifted_chr, fill = list(n = 0))

    # Sort chromosomes numerically where possible
    chr_order <- function(chrs) {
      chr_nums <- str_extract(chrs, "\\d+")
      chr_order <- order(as.numeric(chr_nums), na.last = TRUE)
      chrs[chr_order]
    }

    plot_data <- matrix_data %>%
      mutate(
        icsasg_chr = factor(icsasg_chr, levels = chr_order(unique(icsasg_chr))),
        ssal_lifted_chr = factor(ssal_lifted_chr, levels = chr_order(unique(ssal_lifted_chr)))
      )

    # Highlight selected cell
    selected_data <- NULL
    if (!is.null(data$selected_icsasg_chr) && !is.null(data$selected_ssal_chr)) {
      selected_data <- plot_data %>%
        filter(icsasg_chr == data$selected_icsasg_chr,
               ssal_lifted_chr == data$selected_ssal_chr)
    }

    p <- ggplot(plot_data, aes(x = ssal_lifted_chr, y = icsasg_chr, fill = n)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = ifelse(n > 0, n, "")), size = 3) +
      scale_fill_gradient(low = "white", high = "steelblue",
                          trans = "log1p",
                          breaks = c(0, 1, 10, 100, 1000)) +
      labs(x = "Ssal_v3.1 Chromosome",
           y = "ICSASG_v2 Chromosome",
           fill = "Gene Count",
           title = "Gene Mappings Between Assemblies") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank())

    # Add border around selected cell
    if (!is.null(selected_data) && nrow(selected_data) > 0) {
      p <- p + geom_tile(data = selected_data,
                         aes(x = ssal_lifted_chr, y = icsasg_chr),
                         fill = NA, color = "red", size = 2)
    }

    p
  })

  # Handle matrix clicks
  observeEvent(input$matrix_click, {
    req(filtered_data())

    click <- input$matrix_click

    # Get unique chromosomes and their order
    chr_order <- function(chrs) {
      chr_nums <- str_extract(chrs, "\\d+")
      chr_order <- order(as.numeric(chr_nums), na.last = TRUE)
      chrs[chr_order]
    }

    ssal_chrs <- chr_order(unique(filtered_data()$ssal_lifted_chr))
    icsasg_chrs <- chr_order(unique(filtered_data()$icsasg_chr))

    # Convert click coordinates to chromosome indices
    ssal_idx <- round(click$x)
    icsasg_idx <- round(click$y)

    if (ssal_idx >= 1 && ssal_idx <= length(ssal_chrs) &&
        icsasg_idx >= 1 && icsasg_idx <= length(icsasg_chrs)) {

      data$selected_ssal_chr <- ssal_chrs[ssal_idx]
      data$selected_icsasg_chr <- icsasg_chrs[icsasg_idx]
    }
  })

  # Selected pair text
  output$selected_pair <- renderText({
    if (is.null(data$selected_icsasg_chr) || is.null(data$selected_ssal_chr)) {
      "Click on a cell in the matrix above to view details"
    } else {
      paste0("ICSASG_v2: ", data$selected_icsasg_chr,
             " → Ssal_v3.1: ", data$selected_ssal_chr)
    }
  })

  # Gene count text
  output$gene_count <- renderText({
    if (is.null(data$selected_icsasg_chr) ||
        is.null(data$selected_ssal_chr)) {
      ""
    } else {
      req(filtered_data())
      subset <- filtered_data() %>%
        filter(icsasg_chr == data$selected_icsasg_chr,
               ssal_lifted_chr == data$selected_ssal_chr)

      if (nrow(subset) > 0) {
        class_counts <- subset %>%
          count(class) %>%
          mutate(text = paste0(class, ": ", n)) %>%
          pull(text) %>%
          paste(collapse = ", ")

        paste0(nrow(subset), " genes (", class_counts, ")")
      } else {
        "No genes in this chromosome pair"
      }
    }
  })

  # Dot plot
  output$dotplot <- renderPlot({
    req(filtered_data())
    req(data$selected_icsasg_chr)
    req(data$selected_ssal_chr)

    # Filter data for selected chromosome pair
    plot_data <- filtered_data() %>%
      filter(icsasg_chr == data$selected_icsasg_chr,
             ssal_lifted_chr == data$selected_ssal_chr)

    if (nrow(plot_data) == 0) {
      # Empty plot with message
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No genes map between these chromosomes",
                 size = 6) +
        theme_void()
    } else {
      # Color mapping for classes
      class_colors <- c(
        "Green" = "#2ecc71",
        "Yellow" = "#f39c12",
        "Red" = "#e74c3c",
        "NotMapped" = "#95a5a6",
        "StableIDOnly" = "#3498db"
      )

      ggplot(plot_data, aes(x = icsasg_pos / 1e6, y = ssal_pos / 1e6, color = class)) +
        geom_point(alpha = 0.6, size = 2) +
        scale_color_manual(values = class_colors,
                           breaks = c("Green", "Yellow", "Red", "NotMapped", "StableIDOnly")) +
        labs(x = paste0("ICSASG_v2 ", data$selected_icsasg_chr, " Position (Mb)"),
             y = paste0("Ssal_v3.1 ", data$selected_ssal_chr, " Position (Mb)"),
             color = "Mapping Quality",
             title = paste0("Gene Positions: ", data$selected_icsasg_chr,
                           " → ", data$selected_ssal_chr)) +
        theme_minimal() +
        theme(legend.position = "right")
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
