\docType{data}
\name{airports}
\alias{airports}
\title{European airports statistical data}
\format{An unbalanced panel of 395 Euripean airports in 2008-2012 (1763  observations) on the following 31 variables.

\describe{

\item{ICAO}{Airport ICAO code}
\item{AirportName}{Airport official name}
\item{Country}{Airport's country name}
\item{longitude}{Airport longitude}
\item{latitude}{Airport latitude}
\item{Year}{Observation year}
\item{PAX}{A number of carried passengers}
\item{ATM}{A number of of air transport movements served by an airport}
\item{Cargo}{A total volume of cargo served by an airport}
\item{Population100km}{A number of inhabitants, living in 100 km around an airport}
\item{Population200km}{A number of inhabitants, living in 200 km around an airport}
\item{Island}{1 if an airpiort is located on an island; 0 otherwise}
\item{GDPpc}{Gross domestic product per capita in airport's NUTS3 region}
\item{RevenueTotal}{Airport total revenue}
\item{RevenueAviation}{Airport aviation revenue}
\item{RevenueNonAviation}{Airport non-aviation revenue}
\item{RevenueHandling}{Airport revenue from handling services}
\item{RevenueParking}{Airport revenue from parking services}
\item{EBITDA}{Airport earnings before interest, taxes, depreciation, and amortization}
\item{NetProfit}{Airport net profit}
\item{DA}{Airport depreciation, and amortization}
\item{StaffCount}{A number of staff employed by an airport}
\item{StaffCost}{Airport staff cost}
\item{RunwayCount}{A number of airport runways}
\item{CheckinCount}{A number of airport check-iun facilities}
\item{GateCount}{A number of airport gates}
\item{TerminalCount}{A number of airport terminals}
\item{ParkingSpaces}{A number of airport parking spaces}
\item{RoutesDeparture}{A number of departure routes, served by an airport}
\item{RoutesArrival}{A number of arrival routes, served by an airport}
\item{Routes}{(RoutesDeparture + RoutesArrival)/2}

}}
\source{
\describe{ \item{}{Eurostat (2013). European Statistics
Database, Statistical Office of the European Communities
(Eurostat)} \item{}{Airports' statistical reports(2011)}
\item{}{Open Flights: Airport, airline and route data
http://openflights.org/ (2013-05-31)} \item{}{TDC (2012).
Informe de fiscalizacion de la imputacion por la entidad
"Aeropuertos Espanoles y Navegacion Aerea" (AENA) a cada
uno de los aeropuertos de los ingresos, gastos, e
inversiones correspondientes a la actividad aeroportuaria,
en los ejercicios 2009 y 2010., Tribunal de Cuentas, Spain,
Doc 938.} \item{}{CIESIN, Columbia University. Gridded
Population of the World: Future Estimates (GPWFE). (2005)}
}
}
\description{
The \code{spfrontier} package includes the dataset
\code{airports}, containing information about European
airports infrastructure and traffic statistics in 2011.
}

