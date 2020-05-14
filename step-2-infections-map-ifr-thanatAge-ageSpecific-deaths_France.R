
####
##
## Estimate COVID-19 infections based on mapping Chinese infection fatality rates via thanatological age
## age-specific deaths derived with global deaths distribution by age
## ----take ifr of Salje et al. (2020)
## ----calculate global deaths distribution by age based on data of Dudel et al. (2020) 
##
####

##
### 1. Ungroup reference infection fatality rate with smooth.spline
##


ungrouped_mode_ifr_by_single_age_france_sp <- to_ungroup(to_ungroup=ifr_by_age_france_salje[,2],
                                                        nr_grouped_years=10)

ungrouped_low95_ifr_by_single_age_france_sp <- to_ungroup(to_ungroup=ifr_by_age_france_salje[,3],
                                                         nr_grouped_years=10)

ungrouped_up95_ifr_by_single_age_france_sp <- to_ungroup(to_ungroup=ifr_by_age_france_salje[,4],
                                                        nr_grouped_years=10)

##
### 2. Ungroup remaining life years (ex) & map reference countrie's infection fatality rate via thanatological age
## (this could take a minute...)
##

mapped_mode_ifr_thanatAge <- map_fr_betw_ref_and_coi_thanatAge(deaths=deaths,
						lt_1950_2020=lt_1950_2020,
						ungrouped_cfr_by_single_age_sp=ungrouped_mode_ifr_by_single_age_france_sp,reference="France")

mapped_low95_ifr_thanatAge <- map_fr_betw_ref_and_coi_thanatAge(deaths=deaths,
						lt_1950_2020=lt_1950_2020,
						ungrouped_cfr_by_single_age_sp=ungrouped_low95_ifr_by_single_age_france_sp,reference="France")

mapped_up95_ifr_thanatAge <- map_fr_betw_ref_and_coi_thanatAge(deaths=deaths,
						lt_1950_2020=lt_1950_2020,
						ungrouped_cfr_by_single_age_sp=ungrouped_up95_ifr_by_single_age_france_sp,reference="France")

#
## Adjust for NAs (that may be there in rare cases when values cannot be mapped) 
#

mapped_mode_ifr_thanatAge 

	for(pop in 1:10){
		pos_na <- which(is.na(mapped_mode_ifr_thanatAge[pop,]))
		if(length(pos_na)>0){
		for(pos in 1:length(pos_na)){
			if(pos_na[pos] < 6){
				mapped_mode_ifr_thanatAge[pop,pos_na[pos]] <- min(mapped_mode_ifr_thanatAge[pop,],na.rm=TRUE)
			}
			if(pos_na[pos] >= 6){
				mapped_mode_ifr_thanatAge[pop,pos_na[pos]] <- mapped_mode_ifr_thanatAge[pop,pos_na[pos]-1]
			}
		} ## for pos
		} ## if
	} ## for pop

mapped_mode_ifr_thanatAge 

##

mapped_low95_ifr_thanatAge 

	for(pop in 1:10){
		pos_na <- which(is.na(mapped_low95_ifr_thanatAge[pop,]))
		if(length(pos_na)>0){
		for(pos in 1:length(pos_na)){
			if(pos_na[pos] < 6){
				mapped_low95_ifr_thanatAge[pop,pos_na[pos]] <- min(mapped_low95_ifr_thanatAge[pop,],na.rm=TRUE)
			}
			if(pos_na[pos] >= 6){
				mapped_low95_ifr_thanatAge[pop,pos_na[pos]] <- mapped_low95_ifr_thanatAge[pop,pos_na[pos]-1]
			}
		} ## for pos
		} ## if
	} ## for pop

mapped_low95_ifr_thanatAge 

##

mapped_up95_ifr_thanatAge 

	for(pop in 1:10){
		pos_na <- which(is.na(mapped_up95_ifr_thanatAge[pop,]))
		if(length(pos_na)>0){
		for(pos in 1:length(pos_na)){
			if(pos_na[pos] < 6){
				mapped_up95_ifr_thanatAge[pop,pos_na[pos]] <- min(mapped_up95_ifr_thanatAge[pop,],na.rm=TRUE)
			}
			if(pos_na[pos] >= 6){
				mapped_up95_ifr_thanatAge[pop,pos_na[pos]] <- mapped_up95_ifr_thanatAge[pop,pos_na[pos]-1]
			}
		} ## for pos
		} ## if
	} ## for pop

mapped_up95_ifr_thanatAge 

##
### 3. Calculate lambda_x (age-specific population fraction of people being infected with COVID-19)
### mapping Chinese IFR via thanatological age and age-specific deaths
##

## Get global age distribution of deaths from previous step 1-b

source("Data/global_age_dist_deaths.R")

##

lambda_mode_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths <- get_lambda_thanatAge_globalPattern_ageSpecificDeaths(deaths=deaths,
										global_pattern=global_age_dist_deaths,
										days_observed=str_obs_ahead_2020_selected[1:length(5:ncol(deaths))],
										wom_select_10y=wom_select_10y,
										men_select_10y=men_select_10y,
										cfr_coi_mapped_rc_china_based_on_thanat_x=mapped_mode_ifr_thanatAge)

## lambda_mode_ifr_china_map_thanatAge_globalPattern_ageSpecificDeaths

lambda_low95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths <- get_lambda_thanatAge_globalPattern_ageSpecificDeaths(deaths=deaths,
										global_pattern=global_age_dist_deaths,
										days_observed=str_obs_ahead_2020_selected[1:length(5:ncol(deaths))],
										wom_select_10y=wom_select_10y,
										men_select_10y=men_select_10y,
										cfr_coi_mapped_rc_china_based_on_thanat_x=mapped_low95_ifr_thanatAge)

## lambda_low95_ifr_china_map_thanatAge_globalPattern_ageSpecificDeaths

lambda_up95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths <- get_lambda_thanatAge_globalPattern_ageSpecificDeaths(deaths=deaths,
										global_pattern=global_age_dist_deaths,
										days_observed=str_obs_ahead_2020_selected[1:length(5:ncol(deaths))],
										wom_select_10y=wom_select_10y,
										men_select_10y=men_select_10y,
										cfr_coi_mapped_rc_china_based_on_thanat_x=mapped_up95_ifr_thanatAge)

## lambda_up95_ifr_china_map_thanatAge_globalPattern_ageSpecificDeaths

##
### 4. Calculate different output: sum(lambda_x), I_x, sum(I_x) 
##

output_mode_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths <- get_output_thanatAge_globalPattern_ageSpecificDeaths(lambda_thanat_x_percent_infected=lambda_mode_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths,
									days_observed=str_obs_ahead_2020_selected[1:length(5:ncol(deaths))],
									wom_select_10y=wom_select_10y,
									men_select_10y=men_select_10y) 

## names(output_mode_ifr_china_map_thanatAge_globalPattern_ageSpecificDeaths)
## output_mode_ifr_china_map_thanatAge_globalPattern_ageSpecificDeaths$total_lambda

output_low95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths <- get_output_thanatAge_globalPattern_ageSpecificDeaths(lambda_thanat_x_percent_infected=lambda_low95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths,
									days_observed=str_obs_ahead_2020_selected[1:length(5:ncol(deaths))],
									wom_select_10y=wom_select_10y,
									men_select_10y=men_select_10y) 

## names(output_low95_ifr_china_map_thanatAge_globalPattern_ageSpecificDeaths)
## output_low95_ifr_china_map_thanatAge_globalPattern_ageSpecificDeaths$total_lambda

output_up95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths <- get_output_thanatAge_globalPattern_ageSpecificDeaths(lambda_thanat_x_percent_infected=lambda_up95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths,
									days_observed=str_obs_ahead_2020_selected[1:length(5:ncol(deaths))],
									wom_select_10y=wom_select_10y,
									men_select_10y=men_select_10y) 

## names(output_up95_ifr_china_map_thanatAge_globalPattern_ageSpecificDeaths)
## output_up95_ifr_china_map_thanatAge_globalPattern_ageSpecificDeaths$total_lambda

##
### 5. Visualize total lambda (based on data downloaded on April 18, 2020)
##

##
### modal ifr
##

require(wesanderson)
pal <- c(wes_palette("Darjeeling1"),wes_palette("Darjeeling2"))


##
### 8. Visualize confirmed cases vs estimated infections (based on data downloaded on May 14, 2020)
### all three (mode, high, low) together
### wide format
##

pdf(file="Output/Figure-2_France.pdf", width=18, height=10, family="Times", pointsize=24, onefile=TRUE)

par(fig = c(0,1,0,1), las=1, mai=c(0.6,0.0,1.2,0.2))

 	plot(x=-100,y=-100,xlim=c(0-1400,12500),ylim=c(0,10.5),xlab="",ylab="",cex.main=0.9,
		main="Confirmed cases vs estimated infections, in thousand, as of May 13, 2020",axes=FALSE)

	text(c(10000,11500),c(10.35,10.35),c("Confirmed","Estimated"),pos=3,cex=0.9,col="black",font=2)
	text(c(10000,11500),c(10.0,10.0),c("cases","infections"),pos=3,cex=0.9,col="black",font=2)
	text(-500,10.1,"Quantiles:",pos=4,cex=0.9,col="black",font=2)

	for(pop in 1:length(countries_selected_JHU)){
	
		current_JHU_country <- countries_selected_JHU[pop]
		current_JHU_country_row <- countries_selected_JHU_row[pop]
		
		current_pop_insert <- current_JHU_country
		if(current_pop_insert=="US"){
			current_pop_insert <- "United States of America" 	
		}
		if(current_pop_insert=="Hubei"){
			current_pop_insert <- "China"
		}
		if(current_pop_insert=="Iran"){
			current_pop_insert <- "Iran (Islamic Republic of)"
		}

			## 95% lower IFR:
			current_confirmed <- confirmed[current_JHU_country_row,5:ncol(confirmed)]/1000
			current_infected <- output_low95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths$total_I[pop,]/1000
			current_infected_low <- output_up95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths$total_I[pop,]/1000

			if(current_JHU_country=="Hubei"){
				current_infected <- output_low95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths$total_I[pop,]/1000
			}

			rect(xleft=current_infected_low[length(current_infected_low)],xright=current_infected[length(current_infected)],ybottom=9.25-1*(pop-1),ytop=9.25-1*(pop-1)+0.5,col=grey(0.8))
			text(11000,9.25-1*(pop-1)+0.05,paste(round(current_infected[length(current_infected)],0),"k)",sep=""),pos=4,col=grey(0.8),font=2,cex=0.8)

			if(pop==1){
				text(current_infected[length(current_infected)],9.25-1*(pop-1)+0.45,"0.975",pos=3,col=grey(0.8),font=2,cex=0.9)
			}

			## modal IFR:
			current_confirmed <- confirmed[current_JHU_country_row,5:ncol(confirmed)]/1000
			current_infected <- output_mode_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths$total_I[pop,]/1000
			current_infected_low <- output_up95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths$total_I[pop,]/1000

			if(current_JHU_country=="Hubei"){
				current_infected <- output_mode_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths$total_I[pop,]/1000
			}

			rect(xleft=current_infected_low[length(current_infected_low)],xright=current_infected[length(current_infected)],ybottom=9.25-1*(pop-1),ytop=9.25-1*(pop-1)+0.5,col=grey(0.6))
			text(11000,9.25-1*(pop-1)+0.25+0.2,paste(round(current_infected[length(current_infected)],0),"k",sep=""),pos=4,col=grey(0.6),font=2)

			if(pop==1){
				text(current_infected[length(current_infected)],9.25-1*(pop-1)+0.45,"0.5",pos=3,col=grey(0.6),font=2,cex=0.9)
			}
		
			## 95% upper IFR:
			current_confirmed <- confirmed[current_JHU_country_row,5:ncol(confirmed)]/1000
			current_infected <- output_up95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths$total_I[pop,]/1000

			if(current_JHU_country=="Hubei"){
				current_infected <- output_up95_ifr_france_map_thanatAge_globalPattern_ageSpecificDeaths$total_I[pop,]/1000
			}

			rect(xleft=current_infected[length(current_infected)],xright=current_infected[length(current_infected)],ybottom=9.25-1*(pop-1),ytop=9.25-1*(pop-1)+0.5,col=grey(0.4))
			text(10500,9.25-1*(pop-1)+0.05,paste("(",round(current_infected[length(current_infected)],0),"k, ",sep=""),pos=4,col=grey(0.4),font=2,cex=0.8)

			if(pop==1){
				text(current_infected[length(current_infected)],9.25-1*(pop-1)+0.45,"0.025",pos=3,col=grey(0.4),font=2,cex=0.9)
			}


			lines(x=c(0,6900),y=rep((9.25-1*(pop-1)+0.25),2),col=pal[pop],lty=2,lwd=1)
			rect(xleft=0,xright=confirmed[current_JHU_country_row,ncol(confirmed)]/1000,ybottom=9.15-1*(pop-1),ytop=9.15-1*(pop-1)+0.7,col=pal[pop],border=NA)
			text(10000,9.25-1*(pop-1)+0.25,paste(round(current_confirmed[length(current_infected)],0),"k",sep=""),pos=4,col=pal[pop],font=2)

	}

	axis(side=1,at=seq(0,7000,500),labels=FALSE,lwd=1,pos=0)
	axis(side=1,at=seq(0,10000,1000),labels=TRUE,lwd=3,pos=0)
	axis(side=2,at=seq(0.5,9.5,1),labels=paste(rev(seq(1,10,1)),". ",rev(country_labels),sep=""),lwd=3,pos=0)

dev.off()
