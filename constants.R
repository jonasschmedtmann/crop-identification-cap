#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Master's Thesis - Remote Sensing                      %
# Environmental Engineering - ISA/UL - Lisbon, Portugal %
# (c) 2014 by Jonas Schmedtmann & Manuel Campagnolo     %
#                                                       %
# CONSTANTS SCRIPT                                      %
#                                                       %
# This script contains the constants used in this       %
# project and will be called by the main() function.    %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Caminhos para shapefiles e imagens de satelite
CAMINHO.SHP <<- "../01 Dados/02 Dados espaciais"
CAMINHO.LANDSAT <<- "../01 Dados/03 Imagens satelite/LANDSAT"

#Strings de projeccao de informacao geografica
#Colocar ptLX_e89.gsb na pasta dada por system.file("proj", package = "rgdal")
PROJ4.IGEOE <<- CRS("+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111112 +k=1 +x_0=200000 +y_0=300000 +ellps=intl +units=m +nadgrids=ptLX_e89.gsb +wktext +no_defs")
PROJ4.ETRS <<- CRS("+proj=tmerc +lat_0=39.66825833333333 +lon_0=-8.133108333333334 +k=1 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")
PROJ4.UTM <<- CRS("+proj=utm +zone=29 +ellps=WGS84 +datum=WGS84 +units=m")
ANO <<- c(2005)

PROJ4.IGEOE.LL <<- CRS("+proj=longlat +lat_0=39.66666666666666 +lon_0=-8.131906111111112 +k=1 +x_0=200000 +y_0=300000 +ellps=intl +units=m +nadgrids=ptLX_e89.gsb +wktext +no_defs")
PROJ4.ETRS.LL <<- CRS("+proj=longlat +lat_0=39.66825833333333 +lon_0=-8.133108333333334 +k=1 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")
PROJ4.UTM.LL <<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m")

#Constantes relativas a correccao de imagens de Landsat
BANDAS.LANDSAT <<- c(1,2,3,4,5,7)
ESOLAR.TM5 <<- c(1983,1796,1536,1031,220,83.44)
ESOLAR.ETM7 <<- c(1997,1812,1533,1039,230.8,84.9)
TAMANHO.CELULA <<- 30
#http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20090027884.pdf

#Coordenadas da area de estudo
AREA.X <<- c( 511467,  568131,  573315,  490717,  511467)
AREA.Y <<- c(4361007, 4368546, 4286457, 4286706, 4361007)

#Nomes das culturas
NOMES_CULT <<- c('Permanent grassland',
                 'Forage crops',
                 'Maize',
                 'Rice',
                 'Fallow',
                 'Wheat',
                 'Poor grassland',
                 'Vineyard',
                 'Non used area',
                 'Barley',
                 'Oat',
                 'Olive grove')

COD_CULT <<- c('PGL', 'FOR', 'MAI', 'RIC', 'FAL', 'WHE', 'POG', 'VYA', 'NUA', 'BAR', 'OAT', 'OLI')

ABRV_CULT <<-  c('Perm. Grassland',
                 'Forage crops',
                 'Maize',
                 'Rice',
                 'Fallow',
                 'Wheat',
                 'Poor grassland',
                 'Vineyard',
                 'Non used area',
                 'Barley',
                 'Oat',
                 'Olive grove')

NOMES_CULT_PT <<-  c('Pastagem permanente',
                     'Superfícies forrageiras',
                     'Milho',
                     'Arroz',
                     'Pousio',
                     'Trigo',
                     'Pastagem pobre',
                     'Vinha',
                     'Superfície não usada',
                     'Cevada',
                     'Aveia',
                     'Olival')

ABRV_CULT_PT <<-   c('Pastagem perm.',
                     'Sup. forrageiras',
                     'Milho',
                     'Arroz',
                     'Pousio',
                     'Trigo',
                     'Pastagem pobre',
                     'Vinha',
                     'Sup. n\u{00E3}o usada',
                     'Cevada',
                     'Aveia',
                     'Olival')

COD_CULT_PT <<- c('PAS', 'FOR', 'MIL', 'ARZ', 'POU', 'TRI', 'PAP', 'VIN', 'SNA', 'CEV', 'AVE', 'OLI')



#Coordenadas da area de estudo - PRE DESCOBERTA DE NUVENS
#AREA.X <- c( 490733,  568131,  573315,  490717,  490733)
#AREA.Y <- c(4365213, 4368546, 4286457, 4286706, 4365213)


