#' A function to explore a processed BarSeq run.
#'
#' @param dir The directory containing the summaries.
#' @details This function starts a shiny app to explore the results of a BarSeq run.
#' @keywords BC32 explore results
#' @export
#' @examples
#' exploreBC(getwd())

exploreBC <- function(dir) {
# Load data
	counts <- read.delim('merged_summary_pooled_count.txt', header = T)
	freq <- read.delim('merged_summary_pooled_freq.txt', header = T)

	# Define UI for application
	ui <- shinyUI(fluidPage(

	  # Application title
	  titlePanel('BarSeq Explorer'),

	  # Sidebar with widgets
	  sidebarLayout(
		sidebarPanel(
		  sliderInput(
			'threshold',
			'Threshold for counts',
			1,
			20000,
			1000
		  ),
		  sliderInput(
			'rank',
			'Threshold for rank',
			1,
			100,
			5
		  ),
		  checkboxInput(
			'pool',
			'Pool distribution bargraphs',
			F
		  ),
		  checkboxGroupInput(
			'samples',
			'Select samples',
			names(counts)[2:length(counts)],
			names(counts)[2])
		),

		# Show output plots
		mainPanel(
		  tabsetPanel(
			tabPanel('Distribution', plotOutput('distPlot')),
			tabPanel('TopClones', plotOutput('topPlot')),
			id = 'tabs'
			),
		  actionButton('savePlot', 'Save Plot')
		)
	  )
	))


	# Define server logic
	server <- shinyServer(function(input, output) {

	  output$distPlot <- renderPlot({

		# subset data
		working_data <- data.frame(seq = counts$seq)
		for(i in 1:length(input$samples)) {
		  working_data <- cbind(working_data, counts[, names(counts) == input$samples[i]])
		}
		names(working_data) <- c('seq', input$samples)

		# draw distribution barplot
		working_data <- melt(working_data) %>% filter(value >= input$threshold)
		working_data$value <- as.numeric(working_data$value)
		plot <- ggplot(working_data) +
		  geom_bar(aes(x = seq, y = value, fill = variable), stat = 'identity') +
		  xlab('barcodes') +
		  ylab('counts') +
		  scale_fill_discrete(guide = guide_legend(title = 'Sample')) +
		  theme(axis.text.x = element_blank())
		if(!input$pool) {
		  plot <- plot + facet_grid(variable ~ .)
		} else {}
		plot
	  })

	  output$topPlot <- renderPlot({

		# subset data
		working_data <- data.frame(seq = freq$seq)
		for(i in 1:length(input$samples)) {
		  working_data <- cbind(working_data, freq[, names(freq) == input$samples[i]])
		}
		names(working_data) <- c('seq', input$samples)

		# draw topclones barplot
		topclones <- NULL
		for (i in 1:length(working_data)) {
		  topclones <- merge(topclones, working_data[head(order(working_data[i], decreasing = T), input$rank), c(1, i)], all = T)
		}
		topclones <- topclones[,-2]
		topclones[is.na(topclones)] <- 0
		m <- melt(topclones)
		m <- group_by(m, variable) %>% arrange(desc(value))
		ggplot(m) +
		  geom_bar(aes(x = factor(1), y = value, fill = seq), stat = 'identity') +
		  facet_grid(.~variable) +
		  xlab('') +
		  ylab('Cummulative Frequency') +
		  ylim(c(0, 100)) +
		  theme_bw(base_size = 20) +
		  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
		  scale_fill_discrete(guide = F)

	  })

	  observeEvent(input$savePlot,
		ggsave(
		  paste0(
			paste(input$tabs,
				  paste(input$samples, sep = "", collapse = "_"),
				  input$threshold,
				  input$rank,
				  input$pool,
				  sep = "_"),
			'.png'))
	  )
	})


	# Run the application
	shinyApp(ui = ui, server = server)
}
