#' Spatial Stochastic Frontier 
#' 
#' A set of tools for estimation (MLE) of various spatial specifications of stochastic frontier models
#' @name spfrontier-package
#' @docType package
#' @title Spatial Stochastic Frontier 
#' @author Dmitry Pavlyuk \email{Dmitry.V.Pavlyuk@@gmail.com}
#' @keywords spatial stochastic frontier
NULL

#' European airports statistical data
#' 
#' The \code{spfrontier} package includes the dataset \code{airports}, 
#' containing information about European airports infrastructure and traffic statistics in 2011.
#' 
#' 
#' @name airports
#' @rdname data-airports
#' @docType data
#' @format A dataframe with 113 observations on the following 10 variables.
#'
#' \describe{
#' 
#' \item{ICAO_code}{Airport ICAO code}
#' \item{name}{Airport title}
#' \item{lon}{Airport longitude}
#' \item{lat}{Airport latitude}
#' \item{PAX}{A number of carried passengers}
#' \item{ATM}{A number of of air transport movements served by an airport}
#' \item{freight}{A total volume of freight served by an airport}
#' \item{runways}{A number of airport runways}
#' \item{checkins}{A number of airport check-in facilities}
#' \item{gates}{A number of airport gates}
#' 
#' }
#' @source 
#' \describe{
#' \item{}{Eurostat (2013). European Statistics Database, Statistical Office of the European Communities (Eurostat)}
#' \item{}{Airports' statistical reports(2011)}
#' }
NULL



#' Greece airports statistical data
#' 
#' The \code{spfrontier} package includes the dataset \code{airports}, 
#' containing information about Greece airports infrastructure and traffic statistics in 2011.
#' 
#' 
#' @name airports.greece
#' @rdname data-airports-greece
#' @docType data
#' @format A dataframe with 39 observations on the following 24 variables.
#'
#' \describe{
#' 
#' \item{name}{Airport title}
#' \item{ICAO_code}{Airport ICAO code}
#' \item{lat}{Airport latitude}
#' \item{lon}{Airport longitude}
#' \item{APM_winter}{A number of passengers carried during winter period}
#' \item{APM_summer}{A number of passengers carried during summer period}
#' \item{APM}{A number of passengers carried (winter + summer)}
#' \item{cargo_winter}{A total volume of cargo served by an airport during winter period}
#' \item{cargo_summer}{A total volume of cargo served by an airport during summer period}
#' \item{cargo}{A number volume of cargo served by an airport (winter + summer)}
#' \item{ATM_winter}{A number of air transport movements served by an airport during winter period}
#' \item{ATM_summer}{A number of air transport movements served by an airport during summer period}
#' \item{ATM}{A number of air transport movements served by an airport (winter + summer)}
#' \item{openning_hours_winter}{A total number openning hours during winter period}
#' \item{openning_hours_summer}{A total number openning hours during summer period}
#' \item{openning_hours}{A total number openning hours (winter + summer)}
#' \item{runway_area}{A total area of airport runways}
#' \item{terminal_area}{A total area of airport terminal(s)}
#' \item{parking_area}{A total area of airport parking area}
#' \item{island}{1 if an airpiort is located on an island; 0 otherwise}
#' \item{international}{1 if an airpiort is international; 0 otherwise}
#' \item{mixed_use}{1 if an airpiort is in mixed use; 0 otherwise}
#' \item{WLU}{A total volume of work load units (WLU) served by an airport}
#' \item{NearestCity}{A road network distance between an airport and its nearest city }
#' 
#' }
#' @source 
#' \describe{
#' \item{}{"Airport efficiency and public investment in Greece" (2010) 
#' In Proceeding of the 2010 International Kuhmo-Nectar Conference on Transport Economics, University of Valencia, Spain.}
#' }
NULL



#' Spanish airports statistical data
#' 
#' The \code{spfrontier} package includes the dataset \code{airports}, 
#' containing information about Spanish airports infrastructure and traffic statistics in 2010.
#' 
#' 
#' @name airports.spain
#' @rdname data-airports-spain
#' @docType data
#' @format A dataframe with 43 observations on the following 16 variables.
#'
#' \describe{
#' 
#' \item{ICAO_code}{Airport ICAO code}
#' \item{name}{Airport title}
#' \item{lon}{Airport longitude}
#' \item{lat}{Airport latitude}
#' \item{year}{Year of observation; 2010}
#' \item{APM}{A number of carried passengers}
#' \item{ATM}{A number of of air transport movements served by an airport}
#' \item{freight}{A total volume of freight served by an airport}
#' \item{revenue}{A total airport revenue}
#' \item{DA}{A total deprecation/amortization costs of airport operation}
#' \item{staff_cost}{A total staff cost of airport operation}
#' \item{runways}{A number of airport runways}
#' \item{terminals}{A number of airport terminals}
#' \item{WLU}{A total volume of work load units (WLU) served by an airport}
#' 
#' }
#' @source 
#' \describe{
#' \item{}{TDC (2012). Informe de fiscalizacion de la imputacion por la entidad "Aeropuertos Espanoles y Navegacion  Aerea" (AENA) a cada uno de los aeropuertos de los ingresos, gastos, e inversiones correspondientes a  la actividad aeroportuaria, en los ejercicios 2009 y 2010., Tribunal de Cuentas, Spain, Doc 938.}
#' \item{}{Eurostat (2013). European Statistics Database, Statistical Office of the European Communities (Eurostat)}
#' \item{}{Airports' statistical reports (2011)}
#' }
NULL