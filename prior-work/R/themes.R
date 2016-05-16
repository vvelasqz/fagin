require(grid)
require(scales, quietly=TRUE)

# This is just the theme_bw
rtheme <- function (base_size=12,
                    base_family="",
                    slant.x=FALSE,
                    slant270.x=FALSE,
                    xlab=TRUE,
                    xtitle=TRUE,
                    grid.y=NULL,
                    draw.legend=TRUE,
                    no.title=FALSE,
                    no.legend.title=FALSE
                   ){
    gt <- theme_grey(base_size=base_size, base_family=base_family) %+replace% 
        theme(
            axis.text.x = element_text(size = rel(0.8), colour = 'black'),
            axis.text.y = element_text(size = rel(0.8), colour = 'black'),
            # axis.text.y = element_blank(),
            axis.ticks.x = element_line(colour = "black"), 
            axis.ticks.y = element_blank(),
            legend.key = element_blank(),
            panel.background = element_rect(fill = "white", colour = NA),
            panel.border = element_rect(fill = NA, colour = "grey50"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey80", colour = "grey50"),
            strip.background = element_rect(fill = "grey80", colour = "grey50")
            ) 
    if(slant.x){
        gt <- gt + 
        theme(
            axis.text.x = element_text(angle=330, hjust=0, vjust=1),
            plot.margin=unit(c(0.5,2,0.5,0.5), "cm")
        )
    }
    if(slant270.x){
        gt <- gt + 
        theme(
            axis.text.x = element_text(angle=270, hjust=0, vjust=1)
            # plot.margin=unit(c(0.5,2,0.5,0.5), "cm")
        )
    }
    if(xlab){
        gt <- gt +
        theme(
            axis.text.y = element_text(size = rel(0.8)),
            axis.ticks.y = element_line(colour = "black") 
        )
    }
    if(!xtitle){
        gt <- gt +
        theme(
            axis.title.x = element_blank()
        )
    }
    if(!is.null(grid.y)){
        if(grid.y == 'double'){
            gt <- gt +
            theme(
                panel.grid.major=element_line(colour='grey95', size=c(0.25, 1))
            )
        }
        if(grid.y == 'sparse'){
            gt <- gt +
            theme(
                panel.grid.major=element_line(colour='grey95', size=c(0.25, 0))
            )
        }
        if(grid.y == 'simple'){
            gt <- gt +
            theme(
                panel.grid.major=element_line(colour='grey95', size=0.25)
            )
        }
    }
    if(!draw.legend){
        gt <- gt + theme(legend.position='none')
    }
    if(no.title){
        gt <- gt + theme(plot.title=element_blank())
    }
    if(no.legend.title){
        gt <- gt + theme(legend.title=element_blank())
    }
    return(gt)
}

color.scheme.fill <- function(){
    scale_colour_brewer(palette="Set1")
}

color.scheme.6line <- function(){
    scale_colour_hue(h=c(0,270), l=50, c=200)
}

color.scheme.2line <- function(){
    scale_color_brewer(palette='Set1')
}

logy <- function(){
    scale_y_continuous(
        trans='log2',
        breaks=trans_breaks('log2', function(x) round(2^x))
    )
}

logx <- function(){
    scale_x_continuous(
        trans='log2',
        breaks=trans_breaks('log2', function(x) round(2^x))
    )
}

nox <- function(margin=FALSE){
    out <- theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()
    )
    if(!margin){
        out <- out + 
        theme(
            panel.margin=unit(c(0,0,0,0), 'cm'),
            plot.margin=unit(c(0,0,0,0), 'cm')
        )
    }
    return(out)
}

noy <- function(margin=FALSE){
    out <- theme(
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    )
    return(out)
}

ypercent <- function(){
    scale_y_continuous(labels=percent)
}

color.map <- function(ps=NULL, include.orf=FALSE, n=99, fill=FALSE){

    orf.color <- '#C0C0C0'

    ps <- as.numeric(ps)

    if(setequal(ps, c(19,18,17,16,15,14,13,12,10,9,7,6,3,2,1)) | n==15){
        colmap <- c(
            '#FF00FF', #19
            '#9600E6', #18
            '#0000CC', #17
            '#0000A6', #16
            '#FFFF00', #15
            '#96E600', #14
            '#00CC00', #13
            '#00B300', #12
            '#009900', #10
            '#00FFFF', #9
            '#00D9D9', #7
            '#00BFBF', #6
            '#009999', #3
            '#808080', #2
            '#000000'  #1
        )
    }

    else if(setequal(ps, c(19,18,17,16,12,10,9,7,6,3,2,1)) | n==12){
        colmap <- c(
            '#FF00FF', #19
            '#9600E6', #18
            '#0000CC', #17
            '#0000A6', #16
            '#00FFFF', #12
            '#009999', #10
            '#FFFF00', #9
            '#96E600', #7
            '#00CC00', #6
            '#009900', #3
            '#808080', #2
            '#000000'  #1
        )
    }

    else if(setequal(ps, c(19,16,10,3,2,1)) | n==6){
        colmap <- c(
            '#FF00FF', #19
            '#0000CC', #16
            '#00FFFF', #10
            '#96E600', #3
            '#808080', #2
            '#000000'  #1
        )
    }

    else if(setequal(ps, c(19,1)) | n==2){
        colmap <- c(
            '#FF0000', #19
            '#000000'  #1
        )
    }

    else{
        print("Inappropriate strata")
        return(FALSE)
    }

    if(include.orf){
        colmap <- c(orf.color, colmap)
    }

    if(fill){
        return(scale_fill_manual(values=colmap))
    } else {
        return(scale_colour_manual(values=colmap))
    }
}

color.map2 <- function(ps=NULL, include.orf=FALSE, n=99, fill=FALSE){

    orf.color <- '#E0E0E0'

    ps <- as.numeric(ps)

    if(setequal(ps, c(19,18,17,16,15,14,13,12,10,9,7,6,3,2,1)) | n==15){
        colmap <- c(
            '#FF0000', #19
            '#FF00FF', #18
            '#9600E6', #17
            '#0000A6', #16
            '#0078FF', #15
            '#00C8FF', #14
            '#00FFFF', #13
            '#00D9D9', #12
            '#009999', #10
            '#FF8000', #9
            '#F3F300', #7 
            '#96E600', #6 
            '#009900', #3
            '#808080', #2
            '#000000'  #1
        )
    }

    else if(setequal(ps, c(19,18,17,16,12,10,9,7,6,3,2,1)) | n==12){
        colmap <- c(
            '#FF0000', #19
            '#FF00FF', #18
            '#9600E6', #17
            '#0000FF', #16
            '#00FFFF', #12
            '#00B4E6', #10
            '#FF8000', #9
            '#F3F300', #7 
            '#96E600', #6 
            '#009900', #3
            '#808080', #2
            '#000000'  #1
        )
    }

    else if(setequal(ps, c(19,18,17,16,12,10,9,7,6,3,2,1)) | n==9){
        colmap <- c(
            '#FF00FF', # all alpha
            '#9600E6', # all beta
            '#FF8000', # coiled coil
            '#0000FF', # segregated alpha/beta
            '#00FFFF', # membrane
            '#009999', # multi-domain
            '#FFFF00', # interspersed alpha/beta
            '#96E600', # small
            '#00CC00'  # undefined
        )
    }

    else if(setequal(ps, c(19,16,10,3,2,1)) | n==6){
        colmap <- c(
            '#FF0000', #19
            '#0000FF', #16
            '#00FFFF', #10
            '#96E600', #3
            '#808080', #2
            '#000000'  #1
        )
    }

    else if(setequal(ps, c(19,1)) | n==2){
        colmap <- c(
            '#FF0000', #19
            '#000000'  #1
        )
    }

    else{
        print("Inappropriate strata")
        return(FALSE)
    }

    if(include.orf){
        colmap <- c(orf.color, colmap)
    }

    if(fill){
        return(scale_fill_manual(values=colmap))
    } else {
        return(scale_colour_manual(values=colmap))
    }
}
