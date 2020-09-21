#' Monthly numbers of road traffic accidents with personal injury in BRD  
#'      
#' @format ACCIDENT is a univariate time series of length 528, start  January 1974, frequency = 12   
#' \describe{
#' \item{ACCIDENT}{Monthly  numbers of  road traffic accidents with personal injury}
#' }  
#' @source  < https://www-genesis.destatis.de/genesis//online?operation=table&code=46241-0002&levelindex=0&levelid=1583749114977>
#' @examples 
#' data(ACCIDENT)
#' ## maybe  tsp(ACCIDENT) ; plot(ACCIDENT) 
"ACCIDENT"          
 

#' Alcohol Demand, UK, 1870-1938. 
#' 
#' @format  ALCINCOME is a threevariate time series of length 69 and 3 variables; start 1870, frequency = 1
#' \describe{ 
#' \item{Y}{log consumption per head}
#' \item{Z}{log real income per head} 
#' \item{X}{log real price}
#' }
#' @examples 
#' data(ALCINCOME)
#' ## maybe  tsp(ALCINCOME) ; plot(ALCINCOME) 
#' @source Durbin & Watson (1951)  <https://doi.org/10.1093/biomet/38.1-2.159>
"ALCINCOME"


#' Monthly beer production in Australia: megalitres. Includes ale and stout. Does not include beverages with alcohol percentage less than 1.15.   
#'    
#' @format BEER is a univariate time series of length 476, start January 1956, end Aug 1995, frequency = 12
#' \describe{
#' \item{BEER}{Monthly production of beer in Australia}
#' } 
#' @source R package tsdl  <https://github.com/FinYang/tsdl>
#' @examples 
#' data(BEER)
#' ## maybe  tsp(BEER) ; plot(BEER) 
"BEER" 

#' Weekly number of births in New York     
#'   
#' @format BLACKOUT is a univariate time series of length 313, 1961 -- 1966
#' \describe{
#' \item{BLACKOUT}{Weekly numbers of births  in New York}
#' }
#' @source Izenman, A. J., and Zabell, S. L. (1981)  <https://www.sciencedirect.com/science/article/abs/pii/
#'     0049089X81900181>
#' @examples 
#' data(BLACKOUT)
#' ## maybe  tsp(BLACKOUT) ; plot(BLACKOUT) 
"BLACKOUT"


#' U.S. annual coffee consumption   
#'     
#' @format COFFEE is a univariate time series of length 61; start 1910, frequency = 1 
#' \describe{
#' \item{COFFEE}{annual coffee-consumption USA, logarithmic transformed}
#' }  
#' @source  R package tsdl <https://github.com/FinYang/tsdl>
#' @examples 
#' data(COFFEE)
#' ## maybe  tsp(COFFEE) ; plot(COFFEE) 
"COFFEE" 

    
#' Market value of DAX      
#'     
#' @format DAX is a multivariate time series of length 12180 and 4 variables
#' \describe{
#' \item{DAY}{Day of the week} 
#' \item{MONTH}{Month} 
#' \item{Year}{Year} 
#' \item{DAX30}{Market value} 
#' } 
#' @examples 
#' data(DAX)
#' ## maybe  tsp(DAX) ; plot(DAX) 
"DAX"


#' Incidences of insulin-dependent diabetes mellitus                                                                                                    
#'
#' @format DIABETES is a univariate time series of length 72, start  January 1979, frequency = 12
#' \describe{
#' \item{DIABETES}{Incidences of insulin-dependent diabetes mellitus}
#' }
#' @source Waldhoer, T., Schober, E. and Tuomilehto, J. (1997) <https://www.sciencedirect.com/science/
#'    article/abs/pii/S0895435696003344>
#' @examples 
#' data(DIABETES)
#' ## maybe  tsp(DIABETES) ; plot(DIABETES) 
"DIABETES"


#' Running yield of public bonds in Austria and Germany                                                                                           
#' 
#' @format DOMINANCE is a bivariate time series of length 167:
#' \describe{
#' \item{X}{Interest rate Germany}
#' \item{Y}{Interest rate Austria}
#' }
#' @source Jaenicke, J. and Neck, R. (1996) <https://doi.org/10.17713/ajs.v25i2.555> 
#' @examples 
#' data(DOMINANCE)
#' ## maybe  tsp(DOMINANCE) ; plot(DOMINANCE) 
"DOMINANCE"

                                                                                                               
#' Number of incoming orders for engines   
#'   
#' @format ENGINES is a univariate time series of length 188, start January 1972 frequency = 12
#' \describe{
#' \item{ENGINES}{Incoming orders for enginesn}
#' }   
#' @examples 
#' data(ENGINES)
#' ## maybe  tsp(ENGINES) ; plot(ENGINES)  
"ENGINES" 


#' Portfolio-Insurance-Strategies   
#'  
#' @format FINANCE is a multivariate time series of length 7529:
#' \describe{
#' \item{CPPI}{first Portfolio-Insurance-Strategy}
#' \item{TIPP}{second Portfolio-Insurance-Strategy}
#' \item{StopLoss}{third Portfolio-Insurance-Strategy}
#' \item{SyntheticPut}{fourth Portfolio-Insurance-Strategy}
#' \item{CASH}{money market investment}
#' }
#' @source Dichtl, H. and Drobetz, W. (2011) <doi:10.1016/j.jbankfin.2010.11.012>
#' @examples 
#' data(FINANCE)
#' ## maybe  tsp(FINANCE) ; plot(FINANCE) 
"FINANCE"


#' Germany's gross domestic product adjusted for price changes 
#'   
#' @format GDP is a univariate time series of length 159, start January 1970, frequency = 4
#' \describe{
#' \item{GDP}{Gross domestic product  adjusted for price changes}
#' }
#' @source  <https://www-genesis.destatis.de/genesis//online?operation=table&code=81000-0002&levelindex=0&levelid=1583750132341>
#' @examples 
#' data(GDP)
#' ## maybe  tsp(GDP) ; plot(GDP) 
"GDP"        


#' Germany's gross domestic product, values of Laspeyres index to base 2000 
#'     
#' @format GDPORIG is a univariate time series of length 159, start January 1970, frequency = 4
#' \describe{
#' \item{GDPORIG}{gross domestic product, values of Laspeyres index to the base 2000}
#' }
#' @source <https://www-genesis.destatis.de/genesis//online?operation=table&code=81000-0002&levelindex=0&levelid=1583750132341>
#' @examples 
#' data(GDPORIG)
#' ## maybe  tsp(GDPORIG) ; plot(GDPORIG) 
"GDPORIG"        

                                                                                                                   
#' Cardiac frequency of a patient    
#'    
#' @format HEARTBEAT is a univariate time series of length 30:
#' \describe{
#' \item{HEARTBEAT}{cardiac frequency of a patient}
#' }  
#' @examples 
#' data(HEARTBEAT)
#' ## maybe  tsp(HEARTBEAT) ; plot(HEARTBEAT) 
"HEARTBEAT"

                                                                                                       
#' HSV's position in the first German soccer league       
#'      
#' @format HSV is a univariate time series of length 47:
#' \describe{
#' \item{HSV}{HSV's position in the first German soccer league}
#' }  
#' @source <https://www.transfermarkt.de/hamburger-sv/platzierungen/verein/41>
#' @examples 
#' data(HSV)
#' ## maybe  tsp(HSV) ; plot(HSV) 
"HSV"  

                                                                                                                                      
#' IBM's stock price       
#'     
#' @format IBM is a univariate time series of length 369, start 17 May  1961
#' \describe{
#' \item{IBM}{IBM's daily stock price}
#' }  
#' @source Box, G. E. P.  and Jenkins, G. M. (1970, ISBN: 978-0816210947)  "Time series analysis: forecasting and control"
#' @examples 
#' data(IBM)
#' ## maybe  tsp(IBM) ; plot(IBM) 
"IBM"  


#' Temperature and consumption of ice cream    
#'    
#' @format ICECREAM is a bivariate time series of length 160:  
#' \describe{ 
#' \item{ICE}{consumption of ice cream} 
#' \item{TEMP}{Temperature in Fahrenheit degrees}  
#' }  
#' @source Hand, D. J., et al. (1994, ISBN: 9780412399206)  "A Handbook of Small Data Sets"
#' @examples 
#' data(ICECREAM)
#' ## maybe  tsp(ICECREAM) ; plot(ICECREAM) 
"ICECREAM"


#' Income orders of a company  
#' 
#' @format INORDER is a univariate time series of length 237, start January 1968, frequency =12 
#' \describe{
#' \item{INORDER}{Income orders of a company}
#' }  
#' @examples 
#' data(INORDER)
#' ## maybe  tsp(INORDER) ; plot(INORDER) 
"INORDER"        

                                                                                                                               
#' Daily subsoil water level and precipitation at pilot well Lith     
#'      
#' @format LITH is a bivariate time series of length 1347:
#' \describe{
#' \item{N}{precipitation amount}    
#' \item{G}{water level} 
#' } 
#' @examples 
#' data(LITH)
#' ## maybe  tsp(LITH) ; plot(LITH) 
"LITH"     

 
#' Level of Luteinzing hormone of a cow    
#'     
#' @format LUHORMONE is a bivariate time series of length 29:
#' \describe{
#' \item{T}{Time in minutes}
#' \item{X}{Level of the Luteinzing-hormone}
#' }    
"LUHORMONE" 
  
#' Annual lynx trappings in a region of North-West Canada. Taken from Andrews and Herzberg (1985).      
#'   
#' @format LYNX is a univariate time series of length 114; start 1821 frequency = 1
#' \describe{
#' \item{LYNX}{annual lynx trappings in a region of North-west Canada}
#' } 
#' @source Andrews, D. F. and Herzberg, A. M.  (1985) "Data" <https://www.springer.com/gp/book/9781461295631>
#' @examples 
#' data(LYNX)
#' ## maybe  tsp(LYNX) ; plot(LYNX) 
"LYNX" 


#' Size of populations of lynxes and snow hares 
#'    
#' @format LYNXHARE is a simulated bivariate  time series from a VAR[1]-model  of length 100:  
#' \describe{ 
#' \item{X}{Number of lynxes} 
#' \item{Y}{Number of snow hares}  
#' }  
#' @examples
#' data(LYNXHARE)
"LYNXHARE" 

                                                                                                         
#' Subsoil water level and precipitation at pilot well L921  
#'      
#' @format L921 is a trivariate time series of length 335:
#' \describe{
#' \item{T}{Day}
#' \item{Y}{Water level}
#' \item{Z}{Supplemented water level}
#' } 
#' @examples 
#' data(L921)
#' ## maybe  tsp(L921) ; plot(L921) 
"L921"                     


#' Atmospheric CO2 concentrations (ppmv) derived from in situ air samples collected at Mauna Loa Observatory, Hawaii   
#'    
#' @format MAUNALOA is a univariate time series of length 735; start March 1958, frequency = 12   
#' \describe{
#' \item{MAUNALOA}{CO2-concentration at Mauna Loa}
#' }  
#' @source  Keeling, C. D. ,  Piper, S. C.,  Bacastow,  R. B., Wahlen, M. ,  Whorf, T. P., Heimann,  M.,  and Meijer, H. A.  (2001)  <https://library.ucsd.edu/dc/object/bb3859642r>  
#' @examples 
#' data(MAUNALOA)
#' ## maybe  tsp(MAUNALOA) ; plot(MAUNALOA) 
"MAUNALOA"
 
                                                                                                                         
#' Stock market price of MDAX   
#'       
#' @format  MDAX is a multivariate time series of length 6181 and 4 variables
#' \describe{
#' \item{DAY}{Day of the week} 
#' \item{MONTH}{Month} 
#' \item{YEAR}{Year} 
#' \item{MDAX}{Opening stock market price} 
#' } 
#' @source <https://www.onvista.de/index/MDAX-Index-323547>
#' @examples 
#' data(MDAX)
#' ## maybe  tsp(MDAX) ; plot(MDAX[,3]) 
"MDAX" 

                                                                                                                                                                
#' Melanoma incidence in Connecticut   
#'     
#' @format MELANOM is a  multivariate time series of length 45 and 3 variables
#' \describe{ 
#' \item{POP}{Population} 
#' \item{RATE}{Incidence} 
#' \item{SUN}{Sunspots} 
#' }  
#' @source Andrews, D. F. and Herzberg, A. M.  (1985) "Data" <https://www.springer.com/gp/book/9781461295631>
#' @examples 
#' data(MELANOM)
#' ## maybe  tsp(MELANOM) ; plot(MELANOM[,-1]) 
"MELANOM"  

                                                                                                        
#' Annual trade of muskrat pelts                   
#'     
#' @format MUSKRAT is a univariate time series of length 62; start 1848, frequency = 1 
#' \describe{
#' \item{MUSKRAT}{annual trade of muskrat pelts}
#' }
#' @source <https://archive.uea.ac.uk/~gj/book/data/mink.dat>
#' @examples 
#' data(MUSKRAT)
#' ## maybe  tsp(MUSKRAT) ; plot(MUSKRAT) 
"MUSKRAT" 
 
                                          
#' Amount of an Oxygen isotope      
#'     
#' @format OXYGEN is a  matrix with 164 rows and 2 columns   
#' \describe{ 
#' \item{T}{Time} 
#' \item{D}{DELTA18O}  
#' }  
#' @source Belecher, J., Hampton, J. S., and Tunnicliffe Wilson, T. (1994,  ISSN: 1369-7412)  "Parameterization of Continuous Time Autoregressive Models for Irregularly Sampled Time Series Data"
#' @examples 
#' data(OXYGEN)
#' ## maybe   plot(OXYGEN[,1],OXYGEN[,2],type="l"); rug(OXYGEN[,1])
"OXYGEN" 

                                          
#' Two measurements at a paper machine      
#'     
#' @format PAPER is a bivariate time series of length 160  
#' \describe{ 
#' \item{H}{High} 
#' \item{W}{Weight}  
#' }  
#' @source Janacek, G. J.  & Swift, L. (1993, ISBN: 978-0139184598) "Time Series: Forecasting, Simulation, Applications"
#' @examples 
#' data(PAPER)
#' ## maybe  tsp(PAPER) ; plot(PAPER) 
"PAPER"


#' Monthly prices for pigs  
#'   
#' @format PIGPRICE is a univariate time series of length 240; start January 1894, frequency =12
#' \describe{
#' \item{PIGPRICE}{Monthly prices for pigs}
#' }  
#' @source Hanau, A. (1928)  "Die Prognose der Schweinepreise"
#' @examples 
#' data(PIGPRICE)
#' ## maybe  tsp(PIGPRICE) ; plot(PIGPRICE) 
"PIGPRICE"  

                                                                                                                          
#' Peak power demand in Berlin      
#'     
#' @format PPDEMAND is a univariate time series of length 37; start 1955, frequency = 1
#' \describe{
#' \item{PPDEMAND}{annual peak power demand in Berlin, Megawatt}
#' }  
#' @source Fiedler, H. (1979)  "Verschiedene Verfahren zur Prognose des des Stromspitzenbedarfs in Berlin (West)"
#' @examples 
#' data(PPDEMAND)
#' ## maybe  tsp(PPDEMAND) ; plot(PPDEMAND) 
"PPDEMAND"   


#' Production index of manufacturing industries 
#'     
#' @format PRODINDEX is a univariate time series of length 119:
#' \describe{
#'   \item{PRODINDEX}{Production index of manufacturing industries}
#' } 
#' @source Statistisches Bundesamt (2009)  <https://www-genesis.destatis.de/genesis/online>
#' @examples 
#' data(PRODINDEX)
#' ## maybe  tsp(PRODINDEX) ; plot(PRODINDEX) 
"PRODINDEX"

#' Annual amount of rainfall in Los Angeles    
#'    
#' @format RAINFALL is a univariate time series of length 119; start 1878, frequency = 1
#' \describe{
#' \item{RAINFALL}{Amount of  rainfall in Los Angeles}
#' }  
#' @source LA Times (January 28. 1997)
#' @examples 
#' data(RAINFALL)
#' ## maybe  tsp(RAINFALL) ; plot(RAINFALL) 
"RAINFALL" 

#' Monthly sales of Australian red wine (1000 l)    
#'    
#' @format REDWINE is a univariate time series of length 187; start January 1980, frequency =12
#' \describe{
#' \item{REDWINE}{Monthly sales of Australian red wine }
#' }   
#' @source  R package tsdl <https://github.com/FinYang/tsdl>
#' @examples 
#' data(REDWINE)
#' ## maybe  tsp(REDWINE) ; plot(REDWINE) 
"REDWINE" 
                                                                            
#' CO2-Concentration obtained in Schauinsland, Germany     
#'    
#' @format SCHAUINSLAND is a univariate time series of length 72:
#' \describe{
#' \item{SCHAUINSLAND}{CO2-Concentration obtained in Schauinsland}
#' }  
#' @source <http://cdiac.ornl.gov/trends/co2/uba/uba-sc.html>
#' @examples 
#' data(SCHAUINSLAND)
#' ## maybe  tsp(SCHAUINSLAND) ; plot(SCHAUINSLAND) 
"SCHAUINSLAND" 

#' Annual logging of spruce wood.     
#'    
#' @format SPRUCE is a univariate time series of length 42: 
#' \describe{
#' \item{SPRUCE}{Annual logging of spruce wood}
#' }   
#' @examples 
#' data(SPRUCE)
#' ## maybe  tsp(SPRUCE) ; plot(SPRUCE) 
"SPRUCE"          

                                                                                                  
#' Monthly community taxes in Germany (billions EURO)    
#'    
#' @format TAXES is a univariate time series of length 246; start January 1999, frequency = 12
#' \describe{
#' \item{TAXES}{monthly community taxes in Germany}
#' }   
#' @source   <https://www-genesis.destatis.de/genesis/online?operation=previous&levelindex=1&step=1&titel=Tabellenaufbau&levelid=1583748637039>
#' @examples 
#' data(TAXES)
#' ## maybe  tsp(TAXES) ; plot(TAXES) 
"TAXES"  


#' Mean thickness of annual tree rings 
#'     
#' @format TREERING is a multivariate time series of length 66  with 3 variables:
#' \describe{
#' \item{THICK}{mean thickness of annual tree rings}                                 
#' \item{TEMP}{mean temperature of the year}             
#' \item{RAIN}{amount of rain of the year}
#' } 
#' @source <https://ltrr.arizona.edu/>
#' @examples 
#' data(TREERING)
#' ## maybe  tsp(TREERING) ; plot(TREERING) 
"TREERING"


#' Measurements of physiological tremor    
#'     
#' @format TREMOR is a univariate time series of length 400.
#' \describe{
#' \item{TREMOR}{Tremor}
#' }  
#' @examples 
#' data(TREMOR)
#' ## maybe  tsp(TREMOR) ; plot(TREMOR) 
"TREMOR"                                                                               
    
                                                                              
#' Monthly sales of a company  
#'      
#' @format SALES is a univariate time series of length 77:
#' \describe{
#' \item{y}{monthly  sales of a company}
#' }  
#' @source Newton, H. J. (1988, ISBN: 978-0534091989): "TIMESLAB: A time series analysis laboraty"
#' @examples 
#' data(SALES)
#' ## maybe  tsp(SALES) ; plot(SALES) 
"SALES"   

                                            
#' Population of USA      
#' 
#' @format USAPOP is a univariate time series  of length 39; start 1630, frequency = 0.1 
#' \describe{
#' \item{USAPOP}{Population of USA}
#' }  
#' @source <https://www.worldometers.info/world-population/us-population/>
#' @examples 
#' data(USAPOP)
#' ## maybe  tsp(USAPOP) ; plot(USAPOP) 
"USAPOP" 


#' Concentration of growth hormone of a bull    
#'      
#' @format WHORMONE is a univariate time series of length 97:
#' \describe{
#' \item{WHORMONE}{Concentration of growth hormone of a bull}
#' }  
#' @source Newton, H. J. (1988, ISBN: 978-0534091989): "TIMESLAB: A time series analysis laboraty"
#' @examples 
#' data(WHORMONE)
#' ## maybe  tsp(WHORMONE) ; plot(WHORMONE) 
"WHORMONE"