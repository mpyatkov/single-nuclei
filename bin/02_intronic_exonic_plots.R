#!/usr/bin/env Rscript

library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Intronic/Exonic plot')
    p <- add_argument(p,'--sample_id', default="SAMPLE_ID", help="SAMPLE ID for data")
    p <- add_argument(p,'--input_rds', default="input_rds.rds", help="load precalculated R object")
    p <- add_argument(p, '--number_of_umap_clusters', default=9, help='number of clusters for UMAP plots')
    p <- add_argument(p, '--user_slope', default="auto", help='user defined slope parameter')
    p <- add_argument(p, '--user_intercept', default="auto", help='user defined intercept parameter')
    p <- add_argument(p, '--user_mt_percent', default="10", help='keep only cells less than this MT filtering percent')
    p <- add_argument(p, '--user_population', default="1", help='use only part of the data (1 - above the diagonal, 2- below the diagonal, 3 - both)')
    p <- add_argument(p, '--downstream_gtf', default="genebody", help='Downstream GTF which UMIs will be used')
    return(parse_args(p))
}


argv <- ParseArguments()

# DEBUG <- TRUE
# if (DEBUG){
#   argv$input_rds <- "/projectnb2/wax-dk/max/SCexp/G190_3samples/output/module_1_outputs/rds/G190M2_output_rds.rds"
#   argv$user_slope <- 1
#   argv$user_intercept <- 0
#   argv$user_mt_percent <- 10
#   argv$user_population <- 1
# }


library(gridExtra)
library(grid)
library(patchwork)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(e1071) ## for SVM
library(dplyr)
library(tidyr)
library(stringr)

## "not in" operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)

#BiocManager::install("scDblFinder")# or BiocManager::install("plger/scDblFinder")
library(scDblFinder)

## by default we prefer to use select and filter from dplyr package
select <- dplyr::select
filter <- dplyr::filter

library(purrr)
library(readr)

##### Parameters
# # only for debug
# argv$sample_id <- "SAMPLE_ID"
# argv$user_slope <- "auto"
# argv$user_intercept <- "auto"
# argv$user_mt_gtf <- "genebody"
# argv$user_mt_percent <- "10"
# argv$user_population <- "1" 

###################
#### FUNCTIONS ####
###################

get_gtf_name <- function(param){
    if (param == "with-mono") return("withmono")
    else if (param == "without-mono") return("intronic")
    else if (param == "exonic") return(param)
    else if (param == "genebody") return(param)
    else return(NULL)
}

get_intercepts <- function(df){
    ## each intercept split data on two pieces by percentiles
    ## 33%, mean (not 50% percentile) and 66% percentile
    
    ## mean for exonic_umi and intronic_umi
    mn_e <- mean(df$exonic_umi)
    mn_i <- mean(df$intronic_umi)
    
    ## y = x+intercept. Intercept should be negative for dots which are below diagonal and positive if dots are above the diagonal
    ## this is the main reason why mn_i-mn_e
    real_intercept <- mn_i-mn_e
    
    ## move all points to the origin of the coordinate system 
    df_0 <- data.frame(x=df$exonic_umi-mn_e, y=df$intronic_umi-mn_i) %>% 
        as.matrix()
    
    # rotation matrix for rotating ROW vectors counterclockwise on 45 degree
    mrot_left <- matrix(c(1,-1,1,1), ncol = 2) * matrix(rep(sqrt(2)/2,4), ncol = 2)
    
    # rotate to -45
    df_0_rotated <- df_0 %*% mrot_left
    
    # find percentiles only for the x coordinate because we already rotated dots
    # 33% and 66%
    rot0 <- quantile(df_0_rotated[,1], probs=c(0.33, 0.66))
    
    ## we need to find where is tangent intersect x and y axes
    ## rotate percentiles back to +45 degree by formula 
    ## (0, left) %*% t(mrot_left) -> (new_x, new_y) --> tangent will be 2*new_x  OR the same sqrt(2) * left
    ## (0, right) %*% t(mrot_left) -> ...
    left_percentile <- abs(sqrt(2)*rot0["33%"])
    right_percentile <- abs(sqrt(2)*rot0["66%"])
    
    ## return intercept which allow us to divide the dataset using the diagonal lines with intercepts. 
    ## intercepts split the dataset by percentiles (33%, mean (not 50% percentile), 66%)
    return(round(c(real_intercept+left_percentile, real_intercept, real_intercept-right_percentile), 4))
}


#### PI. FUNCTIONS ####
## convert names: "genebody_05" -> "Genebody - MT 5%"
pretty_names <- function(mt_level) {
    #gtf <- str_extract(mt_level, pattern = "^([:print:]+)(?=_mt)")
    mt_percent <- as.numeric(str_extract(mt_level, pattern = "(?<=_mt)([:print:]+)(?=$)"))

    if (str_starts(mt_level, "genebody")) {
        return(str_interp("Genebody - MT ${mt_percent}%"))
    } else {
        return(str_interp("Intronic + Monoexonic - MT ${mt_percent}%"))
    }
}

## pretty title for each MT plot extracted from panel_01_labels
extract_mito_info <- function(level, labels, full_output = T) {
    ## full_output = T- label with MT cells numbers
    ## full_output = F - label without MT cells numbers, only header
    tmp <- left_join(tibble(mt_level = level), labels, by=c("mt_level")) 
    if (full_output) {
        return(tmp$label)
    } else {
        return(tmp$header)
    }
}
#extract_mito_info("genebody_mt05", panel_01_labels)

## panel_01 -- MT filtering vs GeneBody and Intronic+Monoexonic (separate facet plot for each subset)
plot_one_mt_facet <- function(data, title) {
    ggplot(data, aes(x=exonic_umi, y=intronic_umi, color = as.factor(mt_pass)))+
        geom_point(size=0.1)+
        scale_colour_brewer(labels = c("Mito", "Non-mito"), palette = "Set1")+
        geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
        theme_light()+
        facet_wrap(~mt_pass)+
        guides(color=guide_legend(title="Different levels of MT genes filtering",
                                  override.aes = list(size=5)))+
        theme(legend.position="none",
              legend.text = element_text(size=15),
              legend.title = element_text(size=15),
              strip.background =element_rect(fill="gray90"),
              strip.text = element_text(colour = 'black'),
              plot.title = element_text(size=10))+
        xlab("")+
        ylab("")+
        ggtitle(title)
}

## makes combined plot for Mito and Non-mito cells (mt-pass == {Mito, Non-mito})
plot_one_mt <- function(data) {
    ggplot(data, aes(x=exonic_umi, y=intronic_umi, color = as.factor(mt_pass)))+
        geom_point(size=0.1)+
        scale_colour_brewer(labels = c("Mito", "Non-mito"), palette = "Set1")+
        geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
        theme_light()+
        guides(color=guide_legend(title="Different levels of MT genes filtering",
                                  override.aes = list(size=5)))+
        # theme_light()+
        theme(legend.position="none",
              legend.text = element_text(size=15),
              legend.title = element_text(size=15),
              strip.background =element_rect(fill="gray90"),
              strip.text = element_text(colour = 'black'),
              plot.title = element_text(size=10))+
        ggtitle("Combined mito and non-mito")
}

## combine two previous plots by rows
combined_mt_plot <- function(dt,  title){
    g1 <- plot_one_mt_facet(dt, title)
    g2 <- plot_one_mt(dt)
    cowplot::plot_grid(g1,g2,nrow = 2)+
        theme(panel.border = element_rect(colour = NA, # "black"
                                          fill=NA, 
                                          size=0.1), 
              plot.margin = unit(c(0, 20, 20, 0), "pt"))
}

## input df -->, output list with slope and intercept
svm_slope_and_intercept <- function(df) { ## df --> intronic_umi,exonic_umi,marker
    t_svmfit = svm(marker ~ ., data = df, kernel = "linear", cost = 10, scale = FALSE)
    t_svm_mx <- as.matrix(df %>% select(-marker))
    beta = drop(t(t_svmfit$coefs) %*% t_svm_mx[t_svmfit$index,])
    beta0 = t_svmfit$rho
    list(slope = round(-beta[2]/beta[1], 2),
         intercept = round(beta0/beta[1],2))
}

# aa <- svm_slope_and_intercept(svm_mx)
# p+geom_abline(slope = aa$slope, intercept = aa$intercept)

plot_gtf_based_svm <- function(df, main_title, plot_type="gtf"){
    ## plot_type = {"gtf" | "model"}
    
    ## prepare df for svm
    svm_df <- df %>% 
        select(intronic_umi,exonic_umi, summ2) %>% 
        ## model - all 4 gtf intersection + all dots above the diagonal
        mutate(marker = as.factor(ifelse(summ2 == 4 | intronic_umi > exonic_umi, 1,0))) %>%
        select(-summ2)
    
    if(length(levels(svm_df$marker)) == 1) {
        print(paste0("WARNING: NOT ENOUGH LEVELS FOR SVM"))
        return(NULL)
    }
    
    ## calculate svm
    svm_out <- svm_slope_and_intercept(svm_df)
    
    ## extract information about number of cells
    number_of_cells <- df %>% summarise(n = n()) %>% pull(n)
    
    
    ## calculate number of model cells (Pop.I - 1 and Pop.II - 0)
    num_model_0 <- table(svm_df$marker)[['0']]
    num_model_1 <- table(svm_df$marker)[['1']]
    
    ## adj.population after svm line prediction
    predicted <- svm_df %>%
        mutate(predicted = ifelse(intronic_umi > (svm_out$slope * exonic_umi + svm_out$intercept), "pop1", "pop2" ))
    
    ## return ggplot object up to plot_type
    if (plot_type == "gtf") {
      p <- ggplot(data=df, 
                  aes(x=exonic_umi, 
                      y=intronic_umi, 
                      color=as.factor(summ2)))
    } else {
      p <- ggplot(data=svm_df, 
                  aes(x=exonic_umi, 
                      y=intronic_umi, 
                      color=factor(marker, labels = c("Model pop II", "Model pop I"))))
    }
    
    plot <- if (length(table(predicted$predicted)) == 1) {
      
      print(paste0("WARNING: NOT ENOUGH LEVELS TO PRINT SVM POPULATIONS"))
      ## make title
      t_title_model <- str_interp("${main_title} (No data for SVM)\nNumber of non-mito cells: ${number_of_cells}\nslope=1, intercept=0\nMod.pop I = ${num_model_1}, Mod.pop II = ${num_model_0}\nSVM.pop I = ${num_adj_pop1}, SVM.pop II = ${num_adj_pop2}")
      t_title_gtf <- str_interp("${main_title} (No data for SVM)\nNumber of non-mito cells: ${number_of_cells}\nslope=1, intercept=0")
      
      t_title <- ifelse(plot_type=="gtf", t_title_gtf, t_title_model)
      p +
        geom_point(size=1) +
        geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
        scale_colour_brewer(palette = "Set1")+
        guides(color=guide_legend(title=ifelse(plot_type == "gtf",
                                               "Number of GTFs in point", 
                                               "Model labeling"),
                                  override.aes = list(size=5)))+
        theme_light()+
        theme()+ggtitle(t_title)
    
    } else {
      num_adj_pop1 <- table(predicted$predicted)[['pop1']]
      num_adj_pop2 <- table(predicted$predicted)[['pop2']]
      
      ## make title
      t_title_model <- str_interp("${main_title}\nNumber of non-mito cells: ${number_of_cells}\nslope=${svm_out$slope}, intercept=${svm_out$intercept}\nMod.pop I = ${num_model_1}, Mod.pop II = ${num_model_0}\nAdj.pop I = ${num_adj_pop1}, Adj.pop II = ${num_adj_pop2}")
      t_title_gtf <- str_interp("${main_title}\nNumber of non-mito cells: ${number_of_cells}\nslope=${svm_out$slope}, intercept=${svm_out$intercept}")
      
      t_title <- ifelse(plot_type=="gtf", t_title_gtf, t_title_model)
      
    
      p +
        geom_point(size=1) +
        geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
        geom_abline(intercept=svm_out$intercept, slope =svm_out$slope, color = "blue")+
        scale_colour_brewer(palette = "Set1")+
        guides(color=guide_legend(title=ifelse(plot_type == "gtf",
                                               "Number of GTFs in point", 
                                               "Model labeling"),
                                  override.aes = list(size=5)))+
        theme_light()+
        theme()+ggtitle(t_title)  
    }

    plot
}

plot_one_intercept_slope <- function(df, intercept_column, config = NULL) {
    
    tmp_df <- df
    tmp_table <- df %>% 
        select(!!sym(intercept_column)) %>% 
        table()
    
    pop_1 <- tmp_table['Pop_1']
    pop_2 <- tmp_table['Pop_2']
    
    slope <- as.numeric(str_extract(intercept_column, pattern = "(?<=slope_)([:digit:]+)(?=_intercept)"))
    intercept <- as.numeric(str_extract(intercept_column, pattern = "(?<=intercept_)([:print:]+)(?=$)"))
    
    if (intercept_column == "slope_user_intercept_user") {
        slope = config$slope
        intercept = config$intercept
    }
    
    title <- str_interp("slope = ${slope}, intercept = ${intercept}\nPop. I = ${pop_1}, Pop. II = ${pop_2}")
    
    mcolors <- c("#377eb8", "#e41a1c") # blue, red    
    pop_levels <- c("Population I", "Population II")
    
    ggplot(tmp_df, aes(y = intronic_umi, x = exonic_umi, color = factor(!!sym(intercept_column), levels = c("Pop_1","Pop_2"), lab = c("Population I","Population II"))))+
        geom_point(size = 0.5)+
        geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
        geom_abline(slope=slope, intercept = intercept, color = "red")+
        scale_color_manual(values = mcolors, labels = pop_levels, breaks = pop_levels, limits = pop_levels)+
        # scale_colour_brewer(palette = "Set1", direction = -1)+
        guides(color=guide_legend(title="Populations:",
                                  override.aes = list(size=5)))+
        theme_light() +
        ggtitle(title)
}

#plot_one_intercept_slope(only_filtered, "slope_1_intercept_-0.4386")

## plot 3 slopes plot for given intercept
plot_3slopes <- function(df) {
    
    ## extract string representation of intercepts from columns    
    str_list_intercepts <- names(df) %>% keep(~grepl("slope", .x))
    
    plots <- str_list_intercepts %>% 
        map(function(str_intercept) {
            plot_one_intercept_slope(df, str_intercept)
        })
    
    res <- plots %>% purrr::reduce(`+`) +plot_layout(ncol = 1, nrow = 3) ##, guides = 'collect') & theme(legend.position = 'bottom')
    res
}


#plot_3slopes(only_filtered %>% filter(genebody_mt05 == 1))

## plot panel (1 page) for one mt_level for all intercepts and slopes
plot_all_intercepts <- function(df, mt_level) {
    
    tmp_plots <- list(plot_3slopes(df))
    
    cowplot::plot_grid(purrr::reduce(tmp_plots,`/`) + plot_layout(guides = 'collect') &
                           plot_annotation(title = mt_level)&
                           theme(legend.position = 'bottom'))
}

## add marker trp for specific cell barcodes
dots_in_3GTFs <- function(only_filtered, title) {
    trp_tmp <- only_filtered %>% 
        select(summ, summ2, CB) %>% 
        filter(summ2 == 3) %>% 
        rowwise() %>% 
        mutate(trp = case_when( 
            is.na(summ[[1]]) ~ "Exonic_IntronicMono_GeneBody",
            is.na(summ[[2]]) ~ "Intronic_IntronicMono_GeneBody",
            is.na(summ[[3]]) ~ "Intronic_Exonic_GeneBody",
            TRUE ~ "Intronic_Exonic_IntronicMono")) %>% 
        ungroup()
    
    ## add information about other cell barcodes to plot them as background
    trp <- only_filtered %>% 
        left_join(., trp_tmp %>% select(CB, trp), by="CB") %>% 
        mutate(trp = ifelse(is.na(trp), "Other", trp),
               size = ifelse(trp == "Other", 0.3, 1))
    
    ggplot(trp %>% arrange(size), aes(x = exonic_umi, y = intronic_umi, color = as.factor(trp))) + 
        geom_point(size = 0.5) +
        geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
        theme_light()+
        theme(legend.position = 'none')+
        ggtitle(title)
}

# p <- dots_in_3GTFs(panel_01 %>% filter(mt_pass == "Non-mito", mt_level == "genebody_mt01"), "bla")
# p

panel_03_3gtf_plot <- function(df) {
    
    panel_df <- df %>% 
        filter(mt_pass == "Non-mito") %>% 
        group_by(mt_level) %>%
        nest() %>% 
        mutate(plot = map2(data, mt_level, ~dots_in_3GTFs(., extract_mito_info(mt_level, panel_01_labels, full_output = F)))) %>% 
        arrange(mt_level)
    
    mcolors = brewer.pal(5, name = "Set1")
    mcolors[5] <- "#DADADA" ## gray color for Other cells
    trp_levels = c("Exonic_IntronicMono_GeneBody", 
                   "Intronic_Exonic_GeneBody",
                   "Intronic_Exonic_IntronicMono",
                   "Intronic_IntronicMono_GeneBody", 
                   "Other")
    
    ## patchwork guides='collect' does not work if you have different number of labels in legend
    panel_plot <- purrr::reduce(panel_df$plot,`+`) + 
        plot_layout(ncol = 4, nrow = 2, guides = 'collect') &
        scale_color_manual(values = mcolors, labels = trp_levels, breaks = trp_levels, limits = trp_levels) & 
        guides(color=guide_legend(title="GTF intersections", override.aes = list(size=5))) &
        theme(legend.position = 'bottom', 
              legend.text = element_text(size=11), 
              legend.title = element_text(size=11))
    
    panel_plot
}

## wrapper for all previous functions
plot_01_02_03_panels <- function(panel_01, panel_01_labels, sample_id, downstream_gtf) {
    ## create plots for each mt level separately, each plot represents as ggplot object in 'plot' variable
    panel_01_nested_df <- panel_01 %>% 
        group_by(mt_level) %>%
        nest() %>% 
        mutate(plot = map2(data, mt_level, 
                           ~combined_mt_plot(., extract_mito_info(mt_level, panel_01_labels)))) %>% 
        arrange(mt_level)
    
    # panel_01_nested_df$plot[2]
    panel_01_plot <- purrr::reduce(panel_01_nested_df$plot,`+`)+ plot_layout(ncol = 4, nrow = 2)
    
    # panel_01_plot
    fout <- paste0(sample_id,"_","01_mt_levels_filtering.pdf")
    if (downstream_gtf != "genebody"){
      ggsave(fout, plot = panel_01_plot, width = 18, height = 14)
    }
    
    #### using SVM to separate Non-mito subset of cells as Population I and Population II ########
    ## svm requires that dots already will have a labels
    ## all dots from 4 gtf files and dots where intronic_umi > exonic_umi will be marked as one group
    ## all other dots will be marked as another group
    
    # svm_mx <- panel_01 %>% 
    #   filter(mt_level == "genebody_mt01" & mt_pass == "Non-mito") %>% 
    #   select(intronic_umi,exonic_umi, summ2) %>% 
    #   mutate(marker = as.factor(ifelse(summ2 == 4 | intronic_umi > exonic_umi, 1,0))) %>%
    #   select(-summ2)
    
    ################## panel for model data ####
    panel_02_nested_df_model <- panel_01 %>% 
        filter(mt_pass == "Non-mito") %>% 
        group_by(mt_level) %>%
        nest() %>% 
        mutate(plot = map2(data, 
                           mt_level, 
                           ~plot_gtf_based_svm(., 
                                               extract_mito_info(mt_level, panel_01_labels, full_output = F), 
                                               "model"))) %>% 
        drop_na(plot) %>% 
        arrange(mt_level)
    
    panel_02_plot_model <- purrr::reduce(panel_02_nested_df_model$plot,`+`)+ plot_layout(ncol = 4, nrow = 2, guides = "collect") &
        #guides(color=guide_legend(title="Model labeling", override.aes = list(size=5))) &
        #guides(color=guide_legend(title="Number of GTFs in point", override.aes = list(size=5))) &
        theme(legend.position='bottom',
              legend.text = element_text(size=15),
              legend.title = element_text(size=15)) 
    # panel_02_plot_model
    
    fout <- paste0(argv$sample_id,"_","02_svm_model.pdf")
    #ggsave("panel_02_svm_model.png", plot = panel_02_plot_model, width = 18, height = 14)
    if (downstream_gtf != "genebody"){
      ggsave(fout, plot = panel_02_plot_model, width = 18, height = 14)
    }
    ################## panel for gtf data ####
    panel_02_nested_df_gtf <- panel_01 %>% 
        filter(mt_pass == "Non-mito") %>% 
        group_by(mt_level) %>%
        nest() %>% 
        mutate(plot = map2(data, 
                           mt_level, 
                           ~plot_gtf_based_svm(., 
                                               extract_mito_info(mt_level, panel_01_labels, full_output = F), 
                                               "gtf"))) %>% 
        drop_na(plot) %>% 
        arrange(mt_level)
    
    p2_gtf_colors = brewer.pal(4, name = "Set1")
    p2_gtf_levels = c("1","2","3","4")
    
    panel_02_plot_gtf <- purrr::reduce(panel_02_nested_df_gtf$plot,`+`)+ plot_layout(ncol = 4, nrow = 2, guides = "collect") &
        #guides(color=guide_legend(title="Model labeling", override.aes = list(size=5))) &
        #guides(color=guide_legend(title="Number of GTFs in point", override.aes = list(size=5))) &
        scale_color_manual(values = p2_gtf_colors, labels = p2_gtf_levels, breaks = p2_gtf_levels, limits = p2_gtf_levels) & 
        theme(legend.position='bottom',
              legend.text = element_text(size=15),
              legend.title = element_text(size=15)) 
    #panel_02_plot_gtf
    
    #ggsave("panel_02_svm_gtf.png", plot = panel_02_plot_gtf, width = 18, height = 14)
    fout <- paste0(sample_id,"_","02_svm_gtf.pdf")
    if (downstream_gtf != "genebody"){
      ggsave(fout, plot = panel_02_plot_gtf, width = 18, height = 14)
    }
    #### Plot multiple static slopes and intercepts against different level of contamination filtering ####
    
    # ## specific slopes and intercepts for static slope/intercept plots
    # list_of_intercepts = c(-0.18, -0.15, -0.1, 0)
    
    
    #plot_all_intercepts(only_filtered, "genebody_mt05")
    
    ## plot all (8) panels for each mt_level
    plot_all_mtlevels <- panel_01 %>% 
        filter(mt_pass == "Non-mito") %>% 
        group_by(mt_level) %>%
        nest() %>% 
        mutate(plot = map2(data, mt_level, ~plot_all_intercepts(., extract_mito_info(mt_level, panel_01_labels, full_output = F)))) %>% 
        arrange(mt_level)
    
    ## save all panels as file with multipages
    fout <- paste0(sample_id,"_","02_static_slopes_intercepts.pdf")
    if (downstream_gtf != "genebody"){
    plot_all_mtlevels$plot %>%
        marrangeGrob(nrow = 1, ncol = 1) %>%
        ggsave(filename = fout, width = 14, height = 18)
    }
    #### end
    
    ####
    #### distribution of triples of different GTFs intersections
    ####
    
    
    
    pp <- panel_03_3gtf_plot(panel_01)
    #ggsave("panel_03_3gtf_dots_distribution.png", plot = pp, width = 16, height = 12.40)
    fout <- paste0(sample_id,"_","03_3gtf_dots_distribution.pdf")
    if (downstream_gtf != "genebody"){
      ggsave(fout, plot = pp, width = 16, height = 12.40)
    }
    ### end triple distributions
}

#### P2. FUNCTIONS ####

recalculate_seurat_object <- function(
    seurat_object,
    meta_data,
    config,
    full = T) 
{
        
    ## find optimal number of cluster using binary search for resolution parameter
    binary_search_resolution <- function(expected_clusters, seurat_object) {
        #seurat_object_copy <- seurat_object
        x0 <- 0; x1 <- 1; new_resolution <- 0
        number_of_clusters <- 0
        counter <- 0 
        max_counter <- 15 ## search until max counter otherwise it could be infinite loop for some samples
        while (number_of_clusters != expected_clusters && counter < max_counter)
        {
            new_resolution <- x0+(x1-x0)/2
            
            number_of_clusters <- FindClusters(seurat_object, resolution = new_resolution) %>% 
                .@meta.data %>% 
                select(paste0("SCT_snn_res.",new_resolution)) %>% 
                group_by(paste0("SCT_snn_res.",new_resolution)) %>% 
                distinct() %>%
                summarise(n=n()) %>%
                pull(n)
            
            # print(new_resolution)
            # print(number_of_clusters)
            # print(str_interp("resolution: ${new_resolution} -> number_of_clusters: ${number_of_clusters} x0:${x0} x1:${x1}"))
            
            if (number_of_clusters < expected_clusters) {
                x0 <- new_resolution
            } 
            else {
                x1 <- new_resolution
            }
            
            counter <- counter+1
        }
        
        cat(paste0("==>BINARY SEARCH COUNTER: ",counter,"<=="))
        
        return(new_resolution)
    }
    
    print("==>RUN seurat object<==")
    meta_data <- as.data.frame(meta_data)
    rownames(meta_data) <- meta_data$CB
    new_seurat <- seurat_object[,meta_data$CB]
    ## new_seurat <- subset(seurat_object, cells = meta_data%>%pull(CB))
    print("==>RUN add meta data<==")
    new_seurat <- AddMetaData(new_seurat, meta_data) 
    
    print("==>RUN SCTransform<==")
    
    if (full) {
      new_seurat <- SCTransform(new_seurat, vars.to.regres = config$regression_var)
    }
    
    mito.genes <- grep(pattern = "^(MT|Mt|mt)-", x = rownames(x = new_seurat@assays$SCT@data), value = TRUE)
    var.genes.no <- dplyr::setdiff(new_seurat@assays$SCT@var.features, mito.genes)
    
    print("==>RUN PCA<==")
    n_pcs <- ifelse(nrow(meta_data) < 20, nrow(meta_data)-1, 20)

    n_pcs <- 10
    paste0("NUMBER OF PCS: ", n_pcs) %>% print
    
    if(n_pcs < 20) {
        print(paste0("WARNING: changing numbers of PCs for PCA analysis. Current number of PCs: ", n_pcs))
    }
    
    new_seurat <- RunPCA(new_seurat, npcs = n_pcs, features = var.genes.no)
    
    print("==>RUN UMAP<==")
    n_dims <- ifelse(n_pcs < 10, n_pcs, 10)
    ## https://github.com/satijalab/seurat/issues/4312 (verbose fix)
    new_seurat <- RunUMAP(new_seurat, dims = 1:n_dims,
                          n.neighbors = ifelse(n_dims < 30, n_dims-1, 30L),
                          verbose = ifelse(n_pcs < 50, FALSE, TRUE),
                          min.dist = 0.01)
    
    print("==>RUN FIND NEIGHBOURS<==")
    new_seurat <- FindNeighbors(new_seurat, dims = 1:n_dims)
    
    print("==>RUN OPT RESOLUTION<==")
    opt_resolution <- binary_search_resolution(expected_clusters = config$number_of_umap_clusters, seurat_object = new_seurat)
    print(paste0("Resolution: ",opt_resolution))
    
    print("==>RUN FIND CLUSTERS<==")
    new_seurat <- FindClusters(new_seurat, resolution = opt_resolution)
    
    #DefaultAssay(object = new_seurat) <- "RNA"
    # new_seurat <- new_seurat %>% 
    #     NormalizeData(., assay = "RNA") %>% 
    #     ScaleData(., assay = "RNA")
    #DefaultAssay(object = new_seurat) <- "SCT"
    
    new_seurat
}

seurat_mark_doublets <- function(seurat_obj){
    ## input seurat object without information about singlet/doublet
    ## output seurat object with dblfinder column in meta.data with singlet/doublet
    
    ## do nothing if dblfinder column exists, return object back
    if ("dblfinder" %in% colnames(seurat_obj@meta.data)){
        return(seurat_obj)
    }
    
    if(nrow(seurat_obj@meta.data) < 20) {
        print(paste0("WARNING: only ",nrow(seurat_obj@meta.data)," cellbarcodes available"))
        print("DO NOT CALCULATING DOUBLETS. SET ALL CELLBARCODES AS SINGLETS")
        
        CB <- seurat_obj@meta.data$CB
        dblfinder <- data.frame(CB = CB, dblfinder = "singlet")
        rownames(dblfinder) <- CB
    } else {
        print("==>MARK DOUBLETS<==")
        ## find doublets
        dblfinder <- as.SingleCellExperiment(seurat_obj, assay = "RNA") %>%
            scDblFinder(.) %>%
            as.Seurat(.) %>%
            .$scDblFinder.class %>%
            as.data.frame(.) %>%
            dplyr::select(dblfinder = ".")
    }
    
    new_seurat <- AddMetaData(seurat_obj, dblfinder)
    
    new_seurat
}

seurat_remove_doublets <- function(seurat_obj, config) {
    
    ## mark doublets if it is required    
    new_seurat <- seurat_mark_doublets(seurat_obj)
    
    ## extract doublets/singlets 
    # doublets_singlets <- tibble(CB = names(new_seurat$dblfinder), 
    #                             dblfinder = as.vector(new_seurat$dblfinder))
    print(table(new_seurat@meta.data$dblfinder))
    
    ## assign doublets-singlets to tmp_df to keep only the singlets
    tmp_df <- new_seurat@meta.data %>%
        select(CB, dblfinder) %>% 
        # left_join(., doublets_singlets, by = "CB") %>%
        filter(dblfinder == "singlet")
    
    ## recalculate seurat again with singlets only
    new_seurat <- recalculate_seurat_object(new_seurat,
                                            meta_data =  tmp_df,
                                            config,
                                            full = F)
    new_seurat
}

user_defined_cb_metadata <- function(user_seurat){
    ## output contains CB only for singlets
    
    ## return user clusters, the following datasets based on all user-defined filtering options
    ## if we make left join on whole dataset all NA in corresponding columns will be ("exclude", "doublet")
    meta <- tibble(CB = names(user_seurat$seurat_clusters),
                   user_defined_params_clusters = as.vector(user_seurat$seurat_clusters),
                   user_selected_parameters = "include")
    
    projections <- as_tibble(user_seurat@reductions$umap@cell.embeddings, rownames = "CB") %>% 
        rename("UMAP_1" = "proj_x_user_defined", "UMAP_2" = "proj_y_user_defined")
    
    left_join(meta, projections, by = "CB")
}

export_meta_data <- function(only_filtered, sample_id){
    ## general table with meta data
    meta_data_for_export <- only_filtered %>% 
        select(-summ) %>% 
        rename("intronic_umi" = "log10_intronic_umi",
               "exonic_umi" = "log10_exonic_umi",
               "samples1" = "gtf_intronic",
               "samples2" = "gtf_exonic",
               "samples3" = "gtf_intronic_monoexonic",
               "samples4" = "gtf_genebody",
               "summ2" = "number_of_gtf_intersections",
               "percent.mt_withmono" = "percent.mt_gtf_intronic_monoexonic",
               "percent.mt_genebody" = "percent.mt_gtf_genebody",
               "percent.mt_exonic" = "percent.mt_gtf_exonic",
               "slope_user_intercept_user" = "slope_intercept_user_defined")
    
    ## export all meta.data 
    meta_data_for_export %>% 
        write_csv(paste0(sample_id,"_all_metadata.csv"), col_names = T)
    
    ## export categories (user defined)
    meta_data_for_export %>% 
        select(-log10_intronic_umi, 
               -log10_exonic_umi, 
               -percent.mt_gtf_intronic_monoexonic, 
               -percent.mt_gtf_genebody,
               -all_clusters,
               -proj_x_all,
               -proj_y_all,
               -proj_x_user_defined,
               -proj_y_user_defined) %>% 
        filter(user_selected_parameters == "include") %>% 
        mutate(sample_id = argv$sample_id) %>% 
        write_csv(paste0(sample_id,"_user-defined_cloupe_categories.csv"), col_names = T)
    
    ## export categories (all cellbarcodes defined)
    meta_data_for_export %>% 
        select(-log10_intronic_umi, 
               -log10_exonic_umi, 
               -percent.mt_gtf_intronic_monoexonic, 
               -percent.mt_gtf_genebody,
               -proj_x_all,
               -proj_y_all,
               -proj_x_user_defined,
               -proj_y_user_defined) %>% 
        write_csv(paste0(sample_id,"_all_cloupe_categories.csv"), col_names = T)
    
    ## export projections (user)
    meta_data_for_export %>% 
        filter(user_selected_parameters == "include") %>% 
        select(CB,
               proj_x_user_defined,
               proj_y_user_defined) %>% 
        write_csv(paste0(sample_id,"_user-defined_projections.csv"), col_names = T)
    
    ## export projections (all cellbarcodes)
    meta_data_for_export %>% 
        select(CB,
               proj_x_all,
               proj_y_all) %>% 
        write_csv(paste0(sample_id,"_all_projections.csv"), col_names = T)
    
}

create_umap <- function(seurat_object, lst_entry, config = NULL) {
    
    if (!is.null(config)) {
        ## user defined config
        lst_entry = list(
            slope = config$slope,
            intercept = config$intercept,
            column = "slope_user_intercept_user"
        )
    }
    
    title <- str_interp("slope: ${lst_entry$slope}, intercept: ${lst_entry$intercept}")
    
    if (lst_entry$column == "slope_svm_intercept_svm") {
        title <- str_interp("svm_slope: ${lst_entry$slope}, svm_intercept: ${lst_entry$intercept}")
    }
    
    
    pop1 <- table(seurat_object[[lst_entry$column]])[1]
    pop2 <- table(seurat_object[[lst_entry$column]])[2]
    subtitle <- str_interp("Pop I: ${pop1}, Pop II: ${pop2}")
    
    DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = lst_entry$column)+ggtitle(title) +
        ggtitle(title, subtitle = subtitle) + 
        theme(text=element_text(size=15),
              plot.title = element_text (size=15, hjust = 0),
              plot.subtitle = element_text(size = 15))
}

plot_umap_dotplot <- function(short_seurat) {
  short_seurat <- short_seurat %>% 
    NormalizeData(., assay = "RNA") %>% 
    ScaleData(., assay = "RNA")
  
  DotPlot(short_seurat, features = c(genes_zones, MS_list1, mito), assay = "RNA")+
    RotatedAxis()+
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
          legend.text=element_text(size=7),
          legend.title=element_text(size=7))
}

get_label <- function(df, config) {
    all_nonmito <- df %>% 
        select(starts_with("percent")) %>% 
        mutate(all = n()) %>% 
        filter(percent.mt_genebody < 10) %>% 
        mutate(nonmito = n()) %>% 
        select(all,nonmito) %>% summarize(all = last(all),nonmito = last(nonmito)) %>% as.list(.)
    all <- all_nonmito$all
    nonmito <- all_nonmito$nonmito
    mito <- all - nonmito
    str_interp("${str_to_title(config$short_gtf)} - MT ${config$mt_percent}%, Mito: ${mito}, Non-mito: ${nonmito}")
}

# get mt and all counts per barcode
get_counts_mt_percent <- function(seurat_object, mito, assay = 'RNA') {
    
    # ## per cell barcode
    mt_counts <- as.data.frame(t(seurat_object[[assay]]@counts[mito,])) %>%
        tibble::rownames_to_column(.) %>%
        rowwise(rowname) %>%
        summarise(total = sum(c_across(starts_with("mt")))) %>%
        dplyr::select(CB = rowname, mt_counts = total) %>% 
        ungroup()
    
    # ## per cell barcode
    column <- paste0("nCount_",assay)
    all_counts <- as.data.frame(seurat_object[[column]]) %>%
        as.data.frame(.) %>%
        tibble::rownames_to_column() %>%
        dplyr::select(CB = rowname, all_counts = !!sym(column))
    
    left_join(mt_counts, all_counts, by = "CB")
}

## calculate percent of mt contamination per cluster two different ways
## normalize by rows and by cols (using subset of mito genes)
mt_percent_by_rows_and_columns <- function(seurat_object, mito, regression_var, assay = "RNA") {
    
    ## by rows ([a b, c d] -> mean(a/b + c/d))
    ## return percent for each cluster
    by_rows <- seurat_object@meta.data %>% 
        ## select(CB, seurat_clusters, percent.mt_withmono) %>%
        select(CB, seurat_clusters, !!sym(regression_var)) %>% 
        group_by(seurat_clusters) %>% 
        ## summarise(pt_rows = mean(percent.mt_withmono)) %>% 
        summarise(pt_rows = mean(!!sym(regression_var))) %>% 
        ungroup()
    
    ## by columns ([a b, c d] -> mean((a+c)/(b+d))
    ## return percent for each cluster
    counts_per_barcode <- get_counts_mt_percent(seurat_object = seurat_object, mito = mito, assay = assay)
    
    by_cols <- seurat_object@meta.data %>% 
        select(CB, seurat_clusters) %>% 
        left_join(., counts_per_barcode, by = "CB") %>% 
        mutate(pt = 100*mt_counts/all_counts) %>% 
        group_by(seurat_clusters) %>% 
        summarise(pt_cols = 100*sum(mt_counts)/sum(all_counts),
                  n_cells = n(),
                  MT_UMIs = sum(mt_counts),
                  Total_UMIs = sum(all_counts),
                  `MT_UMIs/cell` = round(MT_UMIs/n_cells, 1),
                  `Total_UMIs/cell` = round(Total_UMIs/n_cells)) %>% 
        ungroup()
    
    left_join(by_rows, by_cols, by = "seurat_clusters") %>% 
        mutate(pt_rows = round(pt_rows,2),
               pt_cols = round(pt_cols,2)) %>% 
        select(cluster = seurat_clusters, pt_rows, pt_cols, n_cells, MT_UMIs, Total_UMIs, `MT_UMIs/cell`, `Total_UMIs/cell`)
}


plot_one_combined_umap <- function(seurat_object, umap_title, config) {
    
    ## extract meta.data for singlets only
    singlet_meta_data <- seurat_object@meta.data %>% 
        select(CB, dblfinder) %>% 
        filter(dblfinder == "singlet")
    
    print(paste0("Number of singlet cells: ", dim(singlet_meta_data)[1]))
    
    ## recalculate clusters after removing doublets
    singlet_seurat_object <- recalculate_seurat_object(seurat_object, 
                                                       meta_data =  singlet_meta_data,
                                                       config,
                                                       full = F)
    ### labels for titles/subtitles
    ## extract number of cells for statistics
    tbl_dblts <- table(seurat_object$dblfinder)
    doublet_num <- tbl_dblts['doublet']
    singlet_num <- tbl_dblts['singlet']
    
    number_of_cells <- length(seurat_object$samples)
    number_of_singlet_cells <- length(singlet_seurat_object$samples)  
    
    ## extract resolution number from "SCT_snn_res.#.###" line for doublet seurat object
    doublets_umap_resolution <- grep("SCT_snn_res.", colnames(seurat_object@meta.data), value = T) %>% 
        gsub("SCT_snn_res.","", .) %>% tail(., n=1)
    
    ## extract resolution number from "SCT_snn_res.#.###" line for singlet seurat object
    singlets_umap_resolution <- grep("SCT_snn_res.", colnames(singlet_seurat_object@meta.data), value = T) %>% 
        gsub("SCT_snn_res.","", .) %>% tail(., n=1)
    
    
    top_umap_subtitle <- str_interp("Number of cells: ${number_of_cells}, umap resolution: ${doublets_umap_resolution}, doublets: ${doublet_num}, singlets: ${singlet_num}")
    
    bottom_umap_subtitle <- str_interp("Number of cells (singlets only): ${number_of_singlet_cells}, umap resolution: ${singlets_umap_resolution}")
    
    ### plots 
    print("Scaling and normalizing data for RNA assay")
    ## test change assay from RNA to SCT
    seurat_object <- NormalizeData(seurat_object, assay = "RNA") %>% 
        ScaleData(., assay = "RNA")
    
    singlet_seurat_object <- NormalizeData(singlet_seurat_object, assay = "RNA") %>% 
        ScaleData(., assay = "RNA")
    
    ## dimplots with doublets, without and with grouping + dotplot
    dimplot_dbl <- DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 7)
    dimplot_dbl_group <- DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 7, group.by = "dblfinder")
    
    ## normalization and scaling only for DotPlot - RNA assay
    seurat_object <- seurat_object %>% 
      NormalizeData(., assay = "RNA") %>%
      ScaleData(., assay = "RNA")
    
    singlet_seurat_object <- singlet_seurat_object %>% 
      NormalizeData(., assay = "RNA") %>%
      ScaleData(., assay = "RNA")
    
    dotplot_dbl <-  DotPlot(seurat_object, features = c(genes_zones, MS_list1, mito), assay = "RNA")+RotatedAxis()
    
    ## dimplot for singlets only + dotplot
    dimplot_sng <- DimPlot(singlet_seurat_object, reduction = "umap", label = TRUE, label.size = 7)
    dotplot_sng <-  DotPlot(singlet_seurat_object, features = c(genes_zones, MS_list1, mito), assay = "RNA")+RotatedAxis()
    
    ## table with mt percents per cluster (assay = RNA- default)
    cluster_percents <- mt_percent_by_rows_and_columns(singlet_seurat_object, mito, regression_var = config$regression_var) 
    cluster_percents_table <- gridExtra::tableGrob(cluster_percents, rows = NULL, theme = ttheme_minimal())
    h <- grobHeight(cluster_percents_table)
    w <- grobWidth(cluster_percents_table)
    title <- textGrob("pt_rows - average mt percent per cluster\npt_cols - 100 * sum MT_UMIs / sum All_UMIs", y=unit(0.5,"npc") + 1.5*h, 
                      vjust=0, gp=gpar(fontsize=20))

    cluster_percents_plot <- gTree(children = gList(title, cluster_percents_table))
        
    top_plot <- cowplot::plot_grid((dimplot_dbl + dimplot_dbl_group) / dotplot_dbl + 
                                       plot_layout(nrow = 2, heights = c(4,1)) + 
                                       plot_annotation(title = umap_title, 
                                                       subtitle = top_umap_subtitle,
                                                       theme = theme(plot.title = element_text(size = 16),
                                                                     plot.subtitle = element_text(size = 12))))
    
    ## layout design for 3 plots (A - UMAP dimplot, B - table, C-dotplot)
    layout_design <- "
  AAAABBB
  AAAABBB
  AAAABBB
  CCCCCCC
  "
    
    bottom_plot <- cowplot::plot_grid(dimplot_sng + cluster_percents_plot + dotplot_sng + 
                                          plot_layout(design = layout_design) + 
                                          plot_annotation(subtitle = bottom_umap_subtitle,
                                                          theme = theme(plot.subtitle = element_text(size = 12)
                                                                        #plot.title = element_text(size = 16),
                                                          )))
    
    cowplot::plot_grid(top_plot/bottom_plot)
    
}

# create plots with user defined parameters (samples_config.csv)
user_defined_plots <- function(df, 
                               downstream_seurat, 
                               config)
{
    
    ## df_user -> slope_user_intercept_user
    ## user_selected_parameters -> include
    
    df_user_mtfiltering <- df %>% filter(!!sym(config$mt_gtf) < config$mt_percent)
    df_user_mtfiltering_population <- df_user_mtfiltering %>% filter(grepl(config$population, slope_user_intercept_user))
    
    ## Prevent execution if not enough cells in dataset
    if(nrow(df_user_mtfiltering_population) < 50) return(list(NULL,NULL))
    
    ## page 1, top plot
    ## real user slope and intercept
    # slope_user_intercept_user <- paste0("slope_",config$slope,"_intercept_",config$intercept)
    # print(slope_user_intercept_user)
    
    top_p1 <- plot_one_intercept_slope(df_user_mtfiltering, "slope_user_intercept_user", config = config)
    
    seurat_user_with_both_populations <- recalculate_seurat_object(downstream_seurat, 
                                                                   meta_data =  df_user_mtfiltering ,
                                                                   USER_DEFINED,
                                                                   full = T)
    
    ## page 1, bottom plot, pop1/pop2 vs umap
    bottom_p1 <- create_umap(seurat_user_with_both_populations, lst_entry = NULL, config)& scale_colour_hue(direction = -1) 
    
    # page 2, whole page plot
    ## TODO: refactor this because somebody would like to change Pop. I and Pop. II to something different
    ## filter by mt percent string samples_config$mt_filtering
    ## filter by population {1,2,3} 1 - I, 2 - II, 3 - I + II
    
    nc_all <-  nrow(df)
    nc_nonmito <- nrow(df_user_mtfiltering) 
    nc_mito <- nc_all - nc_nonmito
    MT_LABEL <- str_interp("sample_id: ${config$sample_id}\nUMIs GTF: ${str_to_title(config$downstream_gtf_name)}, ${str_to_title(config$short_gtf)} - MT ${config$mt_percent}%, Mito: ${nc_mito}, Non-mito: ${nc_nonmito}")
    
    page_title <- str_interp("Only user defined parameters: ${MT_LABEL} mt filtering")
    
    page1 <- cowplot::plot_grid(top_p1/bottom_p1)&plot_annotation(title = page_title, 
                                                                  theme = theme(plot.title = element_text(size = 16),
                                                                                plot.subtitle = element_text(size = 12)))
    ## using only selected population, filtering, singlets
    seurat_only_selected_population <- recalculate_seurat_object(downstream_seurat, 
                                                                 meta_data = df_user_mtfiltering_population,
                                                                 USER_DEFINED,
                                                                 full = T)
    seurat_only_selected_population <- seurat_mark_doublets(seurat_only_selected_population)
    
    page2 <- plot_one_combined_umap(seurat_object = seurat_only_selected_population,
                                    umap_title = page_title,
                                    config = config)
    
    list(patchworkGrob(page1), page2)
}

#### Mapping clusters to Intronic/Exonic plots
combined_clusters_to_scatterplot <- function(seurat_obj_all, seurat_obj_user, df, title, config_default, config_user) {
    
    ## TODO: column for color
    # user_defined_params_clusters
    # all_clusters
    
    title_1 <- str_interp("All CB without filtering. Sample_ID: ${config_default$sample_id}, UMIs GTF: ${str_to_title(config_default$downstream_gtf_name)}\n# of cells ${nrow(df)}")
    umap_all <- DimPlot(seurat_obj_all, reduction = "umap", label = TRUE, label.size = 7)+
        ggtitle(title_1)+
        theme(plot.title = element_text(face = "plain", size = 12))
    
    num_cells <- df %>% filter(user_selected_parameters == "include") %>% nrow()
    pop <- config_user$population
    mtpc <- config_user$mt_percent
    title_2 <- str_interp("CB after applying user defined parameteres. Sample_ID: ${config_user$sample_id}, UMIs GTF: ${str_to_title(config_user$downstream_gtf_name)}\n${str_to_title(config_user$short_gtf)} MT filtering: ${mtpc}%, population: ${pop}, singlets only, # of cells ${num_cells}")
    umap_user <- DimPlot(seurat_obj_user, reduction = "umap", label = TRUE, label.size = 7)+
        ggtitle(title_2)+
        theme(plot.title = element_text(face = "plain", size = 12))
    
    
    left <- ggplot(df, aes(x = exonic_umi, y = intronic_umi, color = all_clusters))+
        geom_point(size = 0.7) + theme_minimal()+
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_abline(slope = config_user$slope, intercept = config_user$intercept, linetype = "solid", color = "blue")+
        guides(color=guide_legend(override.aes = list(size=5), title="All clusters"))+
        ggtitle(str_interp("Clusters without filtering (blue line slope: ${config_user$slope}, intercept: ${config_user$intercept})"))
    
    
    right <- ggplot(df %>% filter(user_selected_parameters == "include") , aes(x = exonic_umi, y = intronic_umi, color = user_defined_params_clusters))+
        geom_point(size = 0.7) + theme_minimal()+
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_abline(slope = config_user$slope, intercept = config_user$intercept, linetype = "solid", color = "blue")+
        guides(color=guide_legend(override.aes = list(size=5), title="User defined"))+
        ggtitle(str_interp("Clusters after applying user defined parameteres (blue line slope: ${config_user$slope}, intercept: ${config_user$intercept})")) 
    
    
    ## normalization and scaling only for DotPlot - RNA assay
    seurat_obj_all <- seurat_obj_all %>% 
      NormalizeData(., assay = "RNA") %>%
      ScaleData(., assay = "RNA")
    
    seurat_obj_user <- seurat_obj_user %>% 
      NormalizeData(., assay = "RNA") %>%
      ScaleData(., assay = "RNA")
    
    all_cb <- DotPlot(seurat_obj_all, 
                      features = c(genes_zones, MS_list1, mito), 
                      assay = "RNA")+
    RotatedAxis()+
      theme(axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
            legend.text=element_text(size=7),
            legend.title=element_text(size=7))
    
    user_cb <- DotPlot(seurat_obj_user, 
                       features = c(genes_zones, MS_list1, mito), 
                       assay = "RNA")+
      RotatedAxis()+
      theme(axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
            legend.text=element_text(size=7),
            legend.title=element_text(size=7))
    
    plot <- cowplot::plot_grid(umap_all + umap_user+left + right  + all_cb + user_cb+
                                   plot_layout(ncol = 2, heights = c(2,2,1)) + 
                                   plot_annotation(title = title, 
                                                   #subtitle = top_umap_subtitle,
                                                   theme = theme(plot.title = element_text(size = 16),
                                                                 plot.subtitle = element_text(size = 12))))

    ## second page
    ## tables with mt percents per cluster
    
    ## Table for all cellbarcodes
    cluster_percents_all <- mt_percent_by_rows_and_columns(seurat_obj_all, mito, regression_var = config_default$regression_var)
    cluster_percents_all_table <- gridExtra::tableGrob(cluster_percents_all, rows = NULL, theme = ttheme_minimal(base_size = 18))
    h <- grobHeight(cluster_percents_all_table)
    w <- grobWidth(cluster_percents_all_table)
    title <- textGrob("Statistics before applying user defined parameters", y=unit(0.5,"npc") + 1.5*h, 
                      vjust=0, gp=gpar(fontsize=20))
    cluster_percents_all_table <- gTree(children = gList(title, cluster_percents_all_table))
    
    ## Table after applying userdefined filtering options
    cluster_percents_user <- mt_percent_by_rows_and_columns(seurat_obj_user, mito, regression_var = config_user$regression_var)
    cluster_percents_user_table <- gridExtra::tableGrob(cluster_percents_user, rows = NULL, theme = ttheme_minimal(base_size = 18))
    h <- grobHeight(cluster_percents_user_table)
    w <- grobWidth(cluster_percents_user_table)
    title <- textGrob("Statistics after applying user defined parameters", y=unit(0.5,"npc") + 1.5*h, 
                      vjust=0, gp=gpar(fontsize=20))
    cluster_percents_user_table <- gTree(children = gList(title, cluster_percents_user_table))
    
    ## TODO: at the moment it is required to create empty ggplot object to plot tables
    tmp_obj <- ggplot() + theme_void()
    
    bottom_plot <- cowplot::plot_grid(tmp_obj+cluster_percents_all_table + cluster_percents_user_table + 
                                          
                                          plot_layout(heights = c(1,10,10)) +
                                          plot_annotation(theme = theme(plot.subtitle = element_text(size = 12))))
    bottom_plot
    
    list(plot, bottom_plot)
}

##############
#### MAIN ####
##############

# List of marker genes for dotplot
genes_zones = c('Cyp2e1', 'Glul', 'Oat', 'Gulo', 'Ass1', 'Hamp', 'Gstp1', 'Ubb',
                'Selenbp2','Sds','Cyp2f2', 'Pck1', 'Hal')

MS_list1 <- c('Il2rb','Cxcr6','Gzma','Csf3r','S100a6','Pkhd1','Sox9','Epcam',
              'Krt7','Krt19','Irf8','Itgax','Clec4f','Csf1r','Jchain','Cd79a',
              'Cd79b','Top2a','Stab2','Kdr','Aqp1','Fcgr2b','Gpr182','Ebf1','Skap1',
              'Ptprc','Ank3','Dcn','Colec11','Ecm1','Alb','Ttr','Apoa1','Serpina1c')

mito <- c("mt-Atp6","mt-Atp8","mt-Co3","mt-Co1","mt-Co2","mt-Nd2","mt-Nd4","mt-Cytb",
          "mt-Nd1","mt-Nd3","mt-Nd4l","mt-Nd5","mt-Nd6")

## load RDS file -> genebody and withmono filtered objects
load(argv$input_rds)
print(argv)

## we need to specify which (genebody or withmono) downstream object
##    and remove the rest

## !!! withmono_raw and genebody_raw in rds contains only union of 4 GTF
## based cell barcodes (not all raw cell barcodes)

downstream_seurat_name <- NULL
downstream_seurat <- NULL

DOWNSTREAM <- argv$downstream_gtf %>% get_gtf_name

if (DOWNSTREAM == "genebody") {
    print("==> FOR DOWNSTREAM ANALYSIS WE WILL USE GENEBODY BASED SEURAT OBJECT <==")
    downstream_seurat_name <- "genebody_raw"
    downstream_seurat <- genebody_raw
} else {
    print("==> FOR DOWNSTREAM ANALYSIS WE WILL USE WITH-MONO BASED SEURAT OBJECT <==")
    downstream_seurat_name <- "withmono_raw"
    downstream_seurat <- withmono_raw
}
rm(genebody_raw,withmono_raw)
gc()

## TODO create downstream_obj
## export
## downstream_seurat_name
## downstream_seurat

# default_mt_gtf <- "withmono"
default_mt_gtf <- DOWNSTREAM
default_mt_percent <- 10
DEFAULT_CONFIG = list(
    downstream_gtf_name = downstream_seurat_name,
    sample_id = argv$sample_id,
    short_gtf = default_mt_gtf,
    mt_gtf = paste0("percent.mt_", default_mt_gtf),
    regression_var = paste0("percent.mt_", default_mt_gtf),
    mt_percent = default_mt_percent,
    ## level for intronic/exonic plot
    mt_gtf_level = paste0(default_mt_gtf,"_mt",default_mt_percent),
    
    population = "Pop_1|Pop_2", # Pop_1, Pop_2, Pop_1|Pop_2
    number_of_umap_clusters = argv$number_of_umap_clusters,
    
    ## TODO: slope can be a vector
    slope = 1,
    intercept = NULL,
    ## gtfs for intronic/exonic plots
    gtfs = c("genebody", "withmono"),
    ## default percent for intronic/exonic plots
    percents = c("05","10","15","20")
)

USER_DEFINED <- list(
    downstream_gtf_name = downstream_seurat_name,
    sample_id = argv$sample_id,
    slope = ifelse(argv$user_slope == "auto", 1, as.numeric(argv$user_slope)),
    intercept = ifelse(argv$user_intercept == "auto", "auto", as.numeric(argv$user_intercept)),
    population = case_when(argv$user_population == "1" ~ "Pop_1",
                           argv$user_population == "2" ~ "Pop_2",
                           TRUE ~ "Pop_1|Pop_2"),
    short_gtf = DOWNSTREAM,
    mt_gtf = paste0("percent.mt_",DOWNSTREAM),
    regression_var = paste0("percent.mt_",DOWNSTREAM),
    mt_percent = as.numeric(argv$user_mt_percent),
    number_of_umap_clusters = argv$number_of_umap_clusters
)

##### End of parameters


# only for debug
## load("./rds_our/G183M1_output_rds_withmono.rds")

## Create column for each (gtf ,percent) combination.
# each column has 'gtf' name and values {0,1} depending on whether 
# the CB passes the filtering threshold  
# https://stackoverflow.com/questions/26003574/use-dynamic-variable-names-in-dplyr
mtfiltering_cols <- cross_df(list(gtf=DEFAULT_CONFIG$gtfs, percent=DEFAULT_CONFIG$percents)) %>% 
    rowwise() %>% 
    mutate(fullname=paste0("percent.mt_",gtf),
           gtf = paste0(gtf,"_mt",percent),
           percent = as.numeric(percent)) %>% 
    pmap_dfc(function(gtf,percent,fullname) { 
        only_filtered %>% 
            select(percent.mt_genebody, percent.mt_withmono) %>% 
            mutate("{gtf}" := ifelse(!!sym(fullname) > percent, 0, 1)) %>% 
            select(!contains("percent"))
    })

## main table which will be used in all downstream calculations
only_filtered <- only_filtered %>% 
    mutate(linear_exonic_umi = exonic_umi,
           linear_intronic_umi = intronic_umi,
           exonic_umi = log10(exonic_umi), 
           intronic_umi = log10(intronic_umi)) %>% 
    bind_cols(., mtfiltering_cols)


## get intercepts for all dots without any filtering 
list_of_intercepts <- get_intercepts(only_filtered)

DEFAULT_CONFIG$intercept <- list_of_intercepts

if (USER_DEFINED$intercept == "auto") {
    ## use 66% percentile as auto intercept for user
    USER_DEFINED$intercept <- list_of_intercepts[3]    
}

rm(list_of_intercepts)

## assing to slope/intercept columns information about population
intercepts_populations <- cross_df(list(slope=DEFAULT_CONFIG$slope, intercept=DEFAULT_CONFIG$intercept)) %>% 
    pmap_dfc(function(slope,intercept) {
        only_filtered %>%
            select(intronic_umi,exonic_umi) %>% 
            mutate("slope_{slope}_intercept_{intercept}" := ifelse(intronic_umi > (exonic_umi * slope + intercept * slope), "Pop_1", "Pop_2")) %>% 
            select(str_interp("slope_${slope}_intercept_${intercept}"))
    }) %>% 
    bind_cols(only_filtered %>% select(CB), .)


only_filtered <- only_filtered %>% left_join(., intercepts_populations, by="CB")
rm(intercepts_populations)

########## PART I. MAKING PLOTS ##########
##### PANEL_01,02,03
## main dataframe which will be used almost for all (excluding umaps) plots on PANEL_01
panel_01 <- only_filtered %>% 
    select(CB, intronic_umi, exonic_umi, summ, summ2, 
           starts_with("genebody"), 
           starts_with("withmono"),
           starts_with("slope")) %>% 
    # transform selected above columns to rows
    pivot_longer( cols = contains("_mt"),
                  names_to = "mt_level",
                  values_to = "mt_pass_value") %>% 
    mutate(mt_pass = ifelse(mt_pass_value == 1, "Non-mito", "Mito"))



## labels for panel 01
panel_01_labels <- panel_01 %>% 
    select(mt_level, mt_pass, mt_pass_value) %>% 
    group_by(mt_level) %>% 
    summarise(non_mito = sum(mt_pass_value),
              mito = n() - non_mito,
              label = paste0("Mito: ", mito, ", Non-mito: ", non_mito)) %>% 
    ungroup() %>% 
    rowwise() %>% 
    mutate(header = pretty_names(mt_level),
           label = str_c(header,", ", label)) %>% 
    select(mt_level, header, label)



plot_01_02_03_panels(panel_01, panel_01_labels, argv$sample_id, DOWNSTREAM)
##### PART II. UMAPs #####

# Elbow plot to select the PC for running umap and finding cluster.
## An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle
## components based on the percentage of variance explained by each one (ElbowPlot() function).
## In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority
## of true signal is captured in the first 10 PCs.
#ElbowPlot(intronic_withmono_raw_v2)


###### Individual umap plots for set of columns which represent ######
# Parameters: slopes: 0.9, 0.95, 1.0
#             intercepts: 0, -0.1, -0.18
# MT_10% - based on intronic+monoexonic                        

## function create subset of the seurat object, all restrictions set in meta_data; 
# number_of_clusters allows user to setup exact number of clusters which will be represented on UMAP plots


## add information about svm and userdefined slope and intercept to only_filtered
## extract slope and intercept for SVM model for intronic+monoexonic with 10% of filtering
svm_params <- panel_01 %>% 
    filter(mt_pass == "Non-mito", mt_level == DEFAULT_CONFIG$mt_gtf_level) %>%
    select(intronic_umi,exonic_umi, summ2) %>% 
    mutate(marker = as.factor(ifelse(summ2 == 4 | intronic_umi > exonic_umi, 1,0))) %>%
    select(-summ2) %>% 
    svm_slope_and_intercept(.)

## add two columns to only filtered
only_filtered <- only_filtered %>% 
    mutate(slope_svm_intercept_svm = ifelse(intronic_umi > (exonic_umi * svm_params$slope + svm_params$intercept), "Pop_1", "Pop_2"),
           slope_user_intercept_user = ifelse(intronic_umi > (exonic_umi * USER_DEFINED$slope + USER_DEFINED$intercept ), "Pop_1", "Pop_2"))

#### SAVE META DATA ####

## applying filtering mt & pop
user_df <- only_filtered %>% 
    filter(!!sym(USER_DEFINED$mt_gtf) < USER_DEFINED$mt_percent) %>% 
    filter(grepl(USER_DEFINED$population, slope_user_intercept_user))

## seurat based on user-defined parameters
user_seurat <- recalculate_seurat_object(downstream_seurat, 
                                         meta_data =  user_df,
                                         USER_DEFINED)
## mark doublets
user_seurat <- seurat_mark_doublets(user_seurat)

## remove doublets
user_seurat_singlets_only <- seurat_remove_doublets(user_seurat, USER_DEFINED)

## only mt filtered, population selected, singlets barcodes 
## returns user_defined_clusters, projections and  "user_selected_parameters columns
user_cb_metadata_df <- user_defined_cb_metadata(user_seurat_singlets_only)

## seurat object with all cellranger filtered cellbarcodes parameters 
all_cb_seurat <- recalculate_seurat_object(downstream_seurat, 
                                           meta_data =  only_filtered, 
                                           DEFAULT_CONFIG)
all_cb_seurat <- seurat_mark_doublets(all_cb_seurat)

## extract only one column with clusters for all cell barcodes from all_cb_seurat
all_cb_metadata_df <- all_cb_seurat@meta.data %>% 
    select(CB, all_clusters = seurat_clusters)

## projections for all cellbarcodes
all_cb_projections <- as_tibble(all_cb_seurat@reductions$umap@cell.embeddings, rownames = "CB") %>% 
    rename("UMAP_1" = "proj_x_all", "UMAP_2" = "proj_y_all")

only_filtered <- only_filtered %>% 
  left_join(., all_cb_metadata_df, by = "CB") %>% 
  left_join(., all_cb_projections, by = "CB") %>% 
  left_join(., user_cb_metadata_df, by = "CB") %>%
  left_join(., user_seurat@meta.data %>% select(CB, dblfinder), by = "CB") %>% 
  replace_na(list(user_selected_parameters = "exclude"))

export_meta_data(only_filtered, argv$sample_id)

############################
### Downstream analysis ####
############################

## df filtered restricted by MT filtering
# df_filtered <- df_filtered %>% 
#     filter(!!sym(SAMPLES_CONFIG$mt_filtering) == 1) ## intronic_mono_mt10

## shortened seurat_object (percent.mt_'gtf' == default_regression_var) 
## NOT FOR USER-DEFINED DOWNSTREAM PLOTS
short_seurat <- recalculate_seurat_object(all_cb_seurat, 
                                          meta_data =  only_filtered %>% 
                                              filter(!!sym(DEFAULT_CONFIG$mt_gtf) < DEFAULT_CONFIG$mt_percent),
                                          DEFAULT_CONFIG,
                                          full = F)
short_seurat <- seurat_mark_doublets(short_seurat)
gc()

## list of "most valuable" slopes and intercepts
slope_intercept_list <- names(only_filtered) %>% 
    keep(~grepl("slope", .x)) %>% 
    discard(~grepl("user",.x)) %>% 
    map(function(r){
        slope <- str_extract(r, pattern = "(?<=slope_)([:print:]+)(?=_intercept)")
        intercept <- str_extract(r, pattern = "(?<=intercept_)([:print:]+)(?=$)")
        
        if (!str_detect(slope,"svm")) {
            slope <- as.numeric(slope)
            intercept <- as.numeric(intercept)
        }
        
        list(real_column = r,
             column = str_replace(r, "-","."),
             slope = slope,
             intercept = intercept)
    })

#create_umap(short_seurat, slope_intercept_list[[3]])

#### create layout for umaps with different slopes and intercepts plots

umap_all_clusters <- DimPlot(short_seurat, reduction = "umap", label = TRUE, label.size = 7)

umap_dot_plot <- plot_umap_dotplot(short_seurat)
gc()

umap_all_clusters_plot <- cowplot::plot_grid(umap_all_clusters,umap_dot_plot,nrow = 2, rel_heights = c(1.95,1.05))

## set layout for 11 plots
layout <- "
AB
CD
EE
EE
"

## create panel 
umap_plots <- slope_intercept_list %>%
    map(function(x) {create_umap(short_seurat, x)}) %>% 
    purrr::reduce(.,`+`) + umap_all_clusters_plot + plot_layout(design = layout, guides = 'collect') &
    scale_colour_hue(direction = -1)&
    guides(color=guide_legend(title="Populations", override.aes = list(size=5))) &
    plot_annotation(title = get_label(only_filtered,DEFAULT_CONFIG)) &
    theme(legend.text = element_text(size=15),
          legend.title = element_text(size=15),
          legend.position = 'top', 
          plot.title = element_text(size = 18))

#umap_plots
fout <- str_glue("{argv$sample_id}_{DOWNSTREAM}_04_umap_plots_populations.pdf")
ggsave(fout, width = 27.14, height = 35)
rm(short_seurat)
##### PANEL_04_2 #####
######## DimPlot and DotPlot for each slope and intercept ####


# t1 <- get_counts_mt_percent(user_seurat_singlets_only, mito = mito, assay = 'RNA')
#t2 <- mt_percent_by_rows_and_columns(user_seurat_singlets_only, mito = mito, assay = 'SCT', regression_var = "percent.mt_genebody")
# t2 <- downstream_seurat %>% NormalizeData(., assay = "RNA") %>% ScaleData(., assay = "RNA")



##### TEST (only specific slope-intercept) ####
# system.time(zz <- plot_one_combined_umap(seurat_object = short_seurat,
#                                          new_meta_data = only_filtered %>% filter(slope_svm_intercept_svm == "Pop_1", withmono_mt10 == 1),
#                                          umap_title = "Intronic+Monoexonic (10% mt filtering)",
#                                          DEFAULT_CONFIG))
# ggsave(plot = zz, filename = "panel_04_umap_plots_combined1.pdf", width = 17.7, height = 22)

#### Make plots for user-defined parameters only ####


#plot_one_intercept_slope(only_filtered, "slope_user_intercept_user", config = USER_DEFINED)

## make 2-pages user defined params plot and save to file
fout <- str_glue("{argv$sample_id}_{DOWNSTREAM}_05_user_defined_only.pdf")
user_defined_plots(df = only_filtered, 
                   downstream_seurat = downstream_seurat, 
                   config = USER_DEFINED) %>% 
    purrr::compact(.) %>% ## do not make plots if they are NULL
    marrangeGrob(nrow = 1, ncol = 1) %>%
    ggsave(filename = fout, width = 17.7, height = 22)

#### Combined UMAP output ####

## system.time(panel_04_umaps2 <- slope_intercept_list %>%
##                 map(function(lst_entry) {
##                   print(paste0("SLOPE_INTERCEPT: ====>",lst_entry$column,"<===="))
##                   gc()
##                   umap_df_filtered <- only_filtered %>%
##                     filter(!!sym(DEFAULT_CONFIG$mt_gtf) < DEFAULT_CONFIG$mt_percent) %>% 
##                     filter(grepl(DEFAULT_CONFIG$population, !!sym(lst_entry$real_column)))
                  
##                   # print(nrow(umap_df_filtered))
##                   ## Prevent execution if not enough cells in dataset
##                   if(nrow(umap_df_filtered) < 50) return(NULL)
##                   tmp_seurat <- recalculate_seurat_object(seurat_object = all_cb_seurat,
##                                                           meta_data = umap_df_filtered %>% select(CB),
##                                                           config = DEFAULT_CONFIG,
##                                                           full = F)
                  
##                   umap_title <- str_interp("${get_label(only_filtered,DEFAULT_CONFIG)} mt filtering) Slope: ${lst_entry$slope}, Intercept: ${lst_entry$intercept}")
##                   plot_one_combined_umap(seurat_object = tmp_seurat,
##                                          umap_title = umap_title,
##                                          config = DEFAULT_CONFIG)
##                 }))

## fout <- paste0(argv$sample_id,"_","04_umap_plots_combined.pdf")
## panel_04_umaps2 %>%
##     purrr::compact(.) %>% 
##     marrangeGrob(nrow = 1, ncol = 1) %>%
##     ggsave(filename = fout, width = 17.7, height = 22)


#### Mapping clusters to Intronic/Exonic plots

rm(user_seurat)

# user_seurat_singlets_only <- user_seurat_singlets_only %>% 
#     NormalizeData(., assay = "RNA") %>% 
#     ScaleData(., assay = "RNA")
# 
# all_cb_seurat <- all_cb_seurat %>% 
#     NormalizeData(., assay = "RNA") %>% 
#     ScaleData(., assay = "RNA")

title <- str_interp("Mapping clusters to intronic/exonic plot. Sample_ID: ${argv$sample_id}")
p1 <- combined_clusters_to_scatterplot(all_cb_seurat, 
                                       user_seurat_singlets_only, 
                                       only_filtered, 
                                       "Mapping clusters intronic/exonic ", 
                                       DEFAULT_CONFIG, 
                                       USER_DEFINED)

fout <- str_glue("{argv$sample_id}_{DOWNSTREAM}_06_clusters_to_scatterplot.pdf")
p1 %>% 
    marrangeGrob(nrow = 1, ncol = 1) %>%
    ggsave(filename = fout, width = 24, height = 18.60)
    # ggsave(filename = fout, width = 17.7, height = 22)


## Temporary stuff
# FeaturePlot(all_cb_seurat, features = c("Alb", "Ttr", "Apoa1","Serpina1c"),cols = c("red","blue"))
# DotPlot(all_cb_seurat, assay = "SCT",features = c(MS_list1, genes_zones, mito))+RotatedAxis()+
#     theme(axis.text.y = element_text(size = 9),
#           axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
#           legend.text=element_text(size=7),
#           legend.title=element_text(size=7))
# 
# FeaturePlot(user_seurat_singlets_only, features = c("Alb", "Ttr", "Apoa1","Serpina1c", "Kdr"))
# 
# DimPlot(user_seurat_singlets_only, reduction = "umap", label = TRUE, label.size = 7, assay = "RNA")
# DotPlot(all_cb_seurat, features = c("Alb", "Ttr", "Apoa1","Serpina1c", "Kdr", "lnc1"), assay = "RNA")
# user_seurat_singlets_only <- NormalizeData(user_seurat_singlets_only, assay = "RNA")
# user_seurat_singlets_only <- ScaleData(user_seurat_singlets_only, assay = "RNA")
# 
# all_cb_seurat <- NormalizeData(all_cb_seurat, assay = "RNA")
# all_cb_seurat <- ScaleData(all_cb_seurat, assay = "RNA")

## TODO: script accidentally creates Rplots.pdf for table with percents, need to fix this latter
## temporary patch: just remove this file 
fn <- "Rplots.pdf"
if (file.exists(fn)) {file.remove(fn)}
