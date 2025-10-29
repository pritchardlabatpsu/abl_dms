library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)

#resistance_data <- read.csv("Highest_resistance_mutation_classification.csv")
resistance_data <- read.csv("Highest_resistance_mutation_classification.csv")
# Known clinical and VUDR variants
known_clinical <- c("T315I", "G250E", "E255K", "M351T", "Y253H", "M244V", "F317L",
                    "F359V", "H396R", "E255V", "E355G", "Q252H", "E459K", "L248V",
                    "Y253F", "D276G", "F359C", "L387M", "F359I", "F486S", "V379I",
                    "L384M", "F311L", "E279K", "H396P", "E355A", "M388L", "F317V",
                    "L387F", "S348L", "D325N", "E255D", "E275K", "F317C", "L248R", "T315V")

vudr_variants <- c("F311I", "Y353H", "K247R", "L298V", "E459G", "S420F", "T277A", "E282Q", "E352G",
                   "E453K", "L273M", "L364I", "R332W", "S417Y", "T315A", "A366V", "A397P", "A407V",
                   "E409K", "E450G", "F382L", "G321E", "M318I", "P441L", "P480L", "R473Q", "T277I",
                   "V289F", "W430L", "Y257C", "Y326H", "A337T", "A365V", "A399T", "A433T", "D276N",
                   "D381Y", "D504E", "E282G", "E292V", "E450A", "E453G", "E453Q", "E466D", "E505K",
                   "E509K", "F317I", "F359A", "F359Y", "G251D", "G254R", "G303E", "G321V", "G426V",
                   "H295Y", "I313L", "I360S", "I418V", "K356R", "M343I", "M351V", "N374S", "P296S",
                   "P309S", "P310L", "P310S", "P310T", "P402S", "P408L", "P441S", "R460H", "S265N",
                   "S385R", "S438C", "S503L", "T267M", "T315N", "Y253C", "Y312C", "Y413C", "A269V",
                   "A288V", "A344T", "A380V", "A397T", "A492D", "C305R", "C305S", "C330G", "C330Y",
                   "C369Y", "D276A", "D325G", "D363G", "D391G", "D391N", "D421N", "D444Y", "D482E",
                   "D482V", "E255L", "E258A", "E258V", "E275Q", "E279A", "E279Y", "E279Z", "E282K",
                   "E282V", "E292D", "E308K", "E316K", "E329D", "E352D", "E352K", "E373G", "E409A",
                   "E409D", "E431G", "E431K", "E450K", "E450Q", "E453A", "E453D", "E453L", "E459Q",
                   "E494G", "E494K", "E499D", "F317R", "F493L", "F497V", "G250R", "G251C", "G251E",
                   "G254E", "G254W", "G303R", "G321R", "G321W", "G383V", "G390E", "G390W", "G398R",
                   "G426K", "H246D", "H246R", "H246Y", "H396Y", "I242M", "I293M", "I313T", "I313V",
                   "I314T", "I360V", "I403T", "I418M", "K247E", "K263N", "K291E", "K291Q", "K357R",
                   "K378R", "K400E", "K404T", "K419E", "L273F", "L284M", "L298R", "L301F", "L324P",
                   "L324Q", "L370R", "L376F", "L411P", "L445M", "M290I", "M290R", "M318T", "M343T",
                   "M351K", "M388I", "M472L", "M496I", "N331S", "N368D", "N368S", "N374Y", "P296L",
                   "P408R", "P439R", "P439S", "P480S", "P484T", "Q252E", "Q252K", "Q252M", "Q252P",
                   "Q252R", "Q300H", "Q346L", "Q447P", "R307L", "R328K", "R328M", "R328W", "R332L",
                   "R332Q", "R362K", "R367Q", "R460C", "R483G", "R483W", "S265R", "S349L", "S438P",
                   "S438Y", "S446P", "S503P", "T243I", "T267A", "T277N", "T315L", "T319I", "T389A",
                   "T394A", "T406I", "T495R", "V256L", "V260A", "V268A", "V268L", "V270M", "V289A",
                   "V299M", "V304G", "V335M", "V338F", "V339A", "V339L", "V371A", "V371I", "V371L",
                   "V377M", "V422I", "W405C", "W423R", "W430C", "W478C", "W478R", "Y253N", "Y264C",
                   "Y264N", "Y320C", "Y440H", "Y449C", "Y469C") 
# Dose color map
dose_colors <- c(
  "400 mg QD" = "#1f78b4",
  "400 mg BID" = "#f46d43",
  "500 mg BID" = "#a50026",
  "Sensitive" = "#f7f7f7"
)

# Full protein × AA grid
full_grid <- expand.grid(
  protein_start = seq(242, 512),
  alt_aa = c("P", "G", "Y", "W", "F", "V", "L", "I", "A", "T",
             "S", "Q", "N", "M", "C", "E", "D", "R", "K", "H")
)

# UI
ui <- fluidPage(
  titlePanel("Imatinib Dose-Dependent Variant Resistance Map"),
  fluidRow(
    column(width = 2,
           wellPanel(
             numericInput("threshold", "Resistance threshold", value = 0.371034195),
             checkboxInput("highlight_clinical", "Highlight Clinical Mutations (green)", TRUE),
             checkboxInput("highlight_vudr", "Highlight VUDRs (magenta)", TRUE),
             checkboxInput("show_legend", "Show Legend", TRUE)
           )
    ),
    column(width = 10,
           plotlyOutput("heatmap_plot", height = "400px", width = "100%")
    )
  )
)

# Server
server <- function(input, output) {
  
  processed_data <- reactive({
    threshold_val <- input$threshold
    
    resistance_data %>%
      mutate(
        Resistant_444nM = rel.viab.444 > threshold_val,
        Resistant_760nM = rel.viab.760 > threshold_val,
        Resistant_916nM = rel.viab.916 > threshold_val,
        Highest_Resistant_Dose = case_when(
          is.na(rel.viab.444) & is.na(rel.viab.760) & is.na(rel.viab.916) ~ NA_character_,
          Resistant_916nM ~ "500 mg BID",
          Resistant_760nM ~ "400 mg BID",
          Resistant_444nM ~ "400 mg QD",
          TRUE ~ "Sensitive"
        ),
        highlight = case_when(
          species %in% known_clinical ~ "Clinical",
          species %in% vudr_variants ~ "Sanger_VUDR",
          TRUE ~ "Normal"
        ),
        tooltip = paste0(
          "Mutation: ", species, "<br>",
          "Highest Resistance dose: ", Highest_Resistant_Dose, "<br>",
          "Viability at 400 mg QD: ", round(rel.viab.444, 3), "<br>",
          "Viability at 400 mg BID: ", round(rel.viab.760, 3), "<br>",
          "Viability at 500 mg BID: ", round(rel.viab.916, 3), "<br>",
          "IC50: ", round(IC50, 3)
        )
      )
  })
  
  merged_data <- reactive({
    full_grid %>%
      left_join(processed_data(), by = c("protein_start", "alt_aa")) %>%
      mutate(fill_value = Highest_Resistant_Dose)
  })
  
  output$heatmap_plot <- renderPlotly({
    merged <- merged_data()
    
    p <- ggplot(merged, aes(x = protein_start, y = alt_aa, fill = fill_value, text = tooltip)) +
      geom_tile(width = 0.95, height = 1, color = NA) +
      
      # Highlight clinical
      { if (input$highlight_clinical)
        geom_tile(data = subset(merged, highlight == "Clinical"),
                  color = "#33a02c", fill = NA, size = 0.3, width = 0.95, height = 1)
      } +
      
      # Highlight VUDR
      { if (input$highlight_vudr)
        geom_tile(data = subset(merged, highlight == "Sanger_VUDR"),
                  color = "#c51b7d", fill = NA, size = 0.3, width = 0.95, height = 1)
      } +
      
      scale_fill_manual(values = dose_colors, na.value = "gray90", name = "Highest Resistance") +
      scale_y_discrete(limits = unique(full_grid$alt_aa)) +
      scale_x_continuous(
        breaks = c(242, 250, 300, 350, 400, 450, 500, 512),
        expand = c(0, 0)
      ) +
      theme_minimal(base_size = 14) +
      labs(x = "Residue on BCR-ABL kinase", y = "Amino Acid Change") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = if (input$show_legend) "bottom" else "none"
      )
    
    ggplotly(p, tooltip = "text")
  })
}

# Run app
shinyApp(ui = ui, server = server)
