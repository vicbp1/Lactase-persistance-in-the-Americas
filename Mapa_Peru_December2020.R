library(raster)
library(rgdal)
library(ggplot2)
library(reshape2)
library(plyr)
library(rgeos)
library(maptools)




nat.earth <- stack("/home/victor/LDGH/Natives_project/Map_natives/NE2_HR_LC_SR/NE2_HR_LC_SR.tif")
ne_boundary <- shapefile("/home/victor/LDGH/Natives_project/Map_natives/ne_50m_admin_0_boundary_lines_land/ne_50m_admin_0_boundary_lines_land.shp")
ne_rivers <- readOGR('/home/victor/LDGH/Natives_project/Map_natives/ne_50m_rivers_lake_centerlines/ne_50m_rivers_lake_centerlines.shp',
                     'ne_50m_rivers_lake_centerlines')
ne_coast <- readOGR('/home/victor/LDGH/Natives_project/Map_natives/ne_50m_coastline/ne_50m_coastline.shp',
                    'ne_50m_coastline')

ne_ocean <- readOGR('/home/victor/LDGH/Natives_project/Map_natives/ne_50m_ocean/ne_50m_ocean.shp',
        'ne_50m_ocean')


quick.subset <- function(x, longlat){
  
  # longlat should be a vector of four values: c(xmin, xmax, ymin, ymax)
  x@data$id <- rownames(x@data)
  
  x.f = fortify(x, region="id")
  x.join = join(x.f, x@data, by="id")
  
  x.subset <- subset(x.join, x.join$long > longlat[1] & x.join$long < longlat[2] &
                       x.join$lat > longlat[3] & x.join$lat < longlat[4])
  
  x.subset
}

domain <- c(-90, -60, -23, 3) ##lon lon lat lat
lakes.subset <- quick.subset(ne_lakes, domain)
boundary.subset <- quick.subset(ne_boundary, domain)
river.subset <- quick.subset(ne_rivers, domain)
coast.subset <- quick.subset(ne_coast, domain)
ocean.subset <- quick.subset(ne_ocean, domain)
nat.crop <- crop(nat.earth, y=extent(domain))

rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                         getValues(nat.crop/255))

rast.table$rgb <- with(rast.table, rgb(NE2_HR_LC_SR.1,
                                       NE2_HR_LC_SR.2,
                                       NE2_HR_LC_SR.3,
                                       1))
  
freq_geo <-  read.table("/home/victor/Projects/Lactose_Intolerance/Mapa/Mapa_Peru_coordinates_frequencies.txt", header = TRUE, sep = "\t")
  
  
  
  ggplot(data = rast.table, aes(x = x, y = y)) +
    geom_tile(fill = rast.table$rgb) +
    geom_path(data=boundary.subset, aes(x = long, y = lat, group = group), size = 0.75,color = 'grey16') +
    geom_point(data=freq_geo, aes(lon, lat, fill=Frequency),color="black", size=7.5,shape=21,stroke=0.8) +
    scale_fill_gradient(low = "white", high ="red" )+ #EDITED
    #theme_void()+
        theme(legend.position = c(-70, -3)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    xlab('') + ylab('')
  
  ggsave("/home/victor/Projects/Lactose_Intolerance/Mapa/Mapa_Peru_Genotypes_Lactase.png", units="cm",width = 19.5, height = 16.9,dpi=300)

  dev.off()
