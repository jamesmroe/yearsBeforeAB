#!/bin/bash 

base=$(echo `dirname \`pwd\``)
echo $base
datdir=$base/data
resdir=$base/results
pltdir=$base/plots
scrdir=$base/scripts
onClust=$2


#---------------------- choose which analysis
analysnum=1			# 1-5
#---------------------- choose analysis step
step=$1					# 1-XX
#-------------------- 


# main analysis opts
updateLCBC=1
updateADNI=1
matchFollow=0				# match groups on age and follow-up data (1 = yes)
procType="LONG"				# longitudinal or cross-sectional freesurfer ("LONG"/"CROSS")
predAB=0					# observed or predicted yearsBeforeAB (1 = predicted)
singleT=0					# all MRI or 1.5T only (1 = 1.5T)


# other set opts
metric="thickness"
savefigs=1
agecut=30
saveres=1
cohortSelect="adnibacslcbc"
yearcuts=$(seq 1 10)
analysisString="tesla"
nTime=2


# set analysis opts by analysis
if [ $analysnum == 2 ]; then
	# sensitivty 1: all MRI + match on follow-up and age
	matchFollow=1
fi

if  [ $analysnum == 3 ]; then
	# sensitivty 2: 1.5T scans only
	singleT=1
fi

if  [ $analysnum == 4 ]; then
	# sensitivty 3: 1.5T scans + match follow-up
	singleT=1; matchFollow=1
fi

if  [ $analysnum == 5 ]; then
	# sensitivty X: cross-sectional proc
	procType="CROSS"
fi

if  [ $analysnum == 6 ]; then
	predAB=1; matchFollow=1
fi

if  [ $analysnum == 7 ]; then
	predAB=1; matchFollow=1; singleT=1
fi


# set analysis name from opts
newstring="reproduce1"
analysname="analysnum${analysnum}matchF${matchFollow}proc${procType}predAB${predAB}singleT${singleT}${newstring}"


# step7 # variable prepSlopes
# step8 # prepslopes compare (main analysis) + resample
# step9 # analyze resampled cortmaps
# step10
# step11


wipeAnalysis=0
if [ $wipeAnalysis == 1 ]; then
	echo "wiping all results for $analysname"
	current_analyses=" \
	cortmaps \
	cortmaps_resample \
	mapPlots \
	prepSlopes \
	prepMods"
	for a in $current_analyses; do
		if [ -d $resdir/$a/$analysname ]; then
			echo "removing $resdir/$a/$analysname"
			rm -r $resdir/$a/$analysname
			echo "removing $base/figures/$analysname"
			rm -r $base/figures/$analysname
		fi
	done
	exit 0
fi



#----------------
# make dirstruct
#----------------
if [ ! -d $datdir ]; then 	
	echo "$(dirname $datdir) does not exist. Quitting";  exit 1
fi
if [ ! -d $resdir ]; then 
	echo "$(dirname $resdir) does not exist. Creating"
	#mkdir $resdir
fi
if [ ! -d $pltdir ]; then 
	echo "$(dirname $pltdir) does not exist. Creating"
	#mkdir $pltdir
fi



if [ $step = 1 ]; then
	echo "----------- preparing slope data -----------"

	for yearsBeforeAB in $yearcuts; do
		
		odir="$resdir/prepSlopes/$analysname"
		ofile="$odir/results.prepSlopes${cohortSelect}${agecut}_${analysisString}yearsBeforeAB${yearsBeforeAB}_${analysname}.Rda"
		# echo $ofile
		
		if [ ! -e "$ofile" ]; then
			
			# set logdir
			logd="logs-01"
			if [ ! -d $logd ]; then 
				mkdir $logd
			fi
			
			if [ $onClust == 1 ]; then
				cmd="sbatch 01-prepSlopesYearsBeforeAB.batch"
			else
				cmd="sh 01-prepSlopesYearsBeforeAB.batch"
			fi
			
			# Run command with args
			$cmd $analysnum $analysname $cohortSelect $yearsBeforeAB $matchFollow $procType $predAB $singleT

		else
			echo "$ofile exists. Skipping."
		fi
	done
fi



if [ $step == 2 ]; then
	echo "----------- running slope comparison (yearsBeforeAB) -----------"
	
	for yearsBeforeAB in $yearcuts; do
		
		odir="$resdir/cortmaps/$analysname/$metric"
		ofile1="$odir/${metric}_cortmapAB${yearsBeforeAB}_nTime${nTime}_${analysname}.csv"
		ofile2="$odir/${metric}_cortmapAB${yearsBeforeAB}_nTime${nTime}_${analysname}_centcor.csv"
		#echo $ofile1; echo $ofile2
		
		if [ ! -e "$ofile1" ] || [ ! -e "$ofile2" ]; then
			logd="logs-02"
			if [ ! -d $logd ]; then
				mkdir $logd
			fi
			
			if [ $onClust == 1 ]; then
				cmd="sbatch 02-prepSlopesCompare.batch"
			else
				cmd="sh 02-prepSlopesCompare.batch"
			fi

			# Run command with args
			$cmd $analysnum $analysname $cohortSelect $yearsBeforeAB $matchFollow $procType $predAB $singleT

		else
			echo "at least one outfile exists: $ofile1. Skipping."
		fi
	done
fi




if [ $step = 3 ]; then
	
	if [ $analysnum -ne 1 ] && [ $analysnum -ne 2 ]; then
		echo "> not running resample. Skip from step 2 to step 5"; exit 0
	fi
	
	echo "----------- running resampled tests (yearsBeforeAB) -----------"
	for yearsBeforeAB in $yearcuts; do

		odir="${base}/results/cortmaps_resample/${analysname}/${metric}/yearsBeforeAB${yearsBeforeAB}nTime${nTime}"
		if [ ! -d $odir ]; then 
			echo "$odir does not exist. Cannot run resample analysis. Quitting";  exit 1
		fi
		logd="logs-03"
		if [ ! -d $logd ]; then
			mkdir $logd
		fi
		# nresamples * 10 yearsBeforeAB
		subset_size=20
		total=$((1000))
		splits=$(((total + subset_size - 1) / subset_size))
		
		for i in $(seq 1 $splits); do
			start=$(( (i * subset_size) - subset_size + 1 ))
			end=$(( start + subset_size - 1 ))
			echo "AB${yearsBeforeAB} Start: $start, End: $end"

			ofile1="$odir/results_cortmap_yearsBeforeAB${yearsBeforeAB}nTime${nTime}-`seq -f "%04g" ${end} ${end}`.csv"
			ofile2="$odir/results_cortmapcentcor_yearsBeforeAB${yearsBeforeAB}nTime${nTime}-`seq -f "%04g" ${end} ${end}`.csv"
			# echo $ofile1; echo $ofile2; exit 0

			if [[ ! -e $ofile1 ]] || [[ ! -e $ofile2 ]]; then
				logd="logs-03"
				if [ ! -d $logd ]; then
					mkdir $logd
				fi
				
				if [ $onClust == 1 ]; then
					cmd="sbatch 03-cortmaps_resample_nTime.batch"
				else
					cmd="sh 03-cortmaps_resample_nTime.batch"
				fi
				echo -n "Submitting job #$i... "
				$cmd $start $end $yearsBeforeAB $analysname $metric $nTime
			else
				echo "at least one outfile exists: $ofile1. Skipping."
			fi
		done
	done
fi


if [ $step = 4 ]; then

	if [ $analysnum -ne 1 ] && [ $analysnum -ne 2 ]; then
		echo "> not running resample. Skip from step 2 to step 5"; exit 0
	fi

	echo "----------- creating summary cortical maps from resampled tests (yearsBeforeAB) -----------"
	
	#NB! loops yearsbefore AB in script
	# yearsBeforeAB=1 #example for ofile
	for yearsBeforeAB in $yearcuts; do
		logd="logs-04"
		if [ ! -d $logd ]; then
			mkdir $logd
		fi
		odir="$resdir/cortmaps/$analysname/$metric"
		ofile="$odir/cortmap_${metric}_resamplesummary_AB${yearsBeforeAB}nTime${nTime}_${analysname}.csv"
		
		if [ ! -e "$ofile" ]; then
			if [ $onClust == 1 ]; then
				cmd="sbatch 04-cortmaps_summarize.batch"
			else
				cmd="sh 04-cortmaps_summarize.batch"
			fi
			$cmd $analysname $yearsBeforeAB $metric $nTime
		else
			echo "$ofile exists. Skipping."
		fi
	done
fi


# plotting procedure only

# if [ $step = 5 ]; then
# 	echo "----------- making ggseg maps (yearsBeforeAB) -----------"
	
# 	if [ $analysnum = 1 ] || [ $analysnum = 2 ]; then
# 		#run also for resmple summary maps
# 		logical="TRUE FALSE"
# 	else
# 		logical="FALSE"
# 	fi
	
# 	if [ $newstring = "longcombat" ] || [ $newstring = "nosingleabminus" ]; then
# 		logical="FALSE"
# 	fi

# 	odir="$resdir/mapPlots/$analysname/$metric"
# 	# echo $odir; exit 0
# 	if [ -d "$odir" ]; then
# 		echo "> removing *pdf and *png already there"
# 		echo "> current fix for making sure fig concat works - rerun map creation when adding new"
# 		echo "$odir"
# 		rm -r $odir
# 	fi
# 	logd="logs-05"
# 	if [ ! -d $logd ]; then
# 		mkdir $logd
# 	fi
# 	empiricalMap="FALSE"
# 	for yearsBeforeAB in $yearcuts; do
# 		for summaryMap in $logical; do
			
# 			if [ $summaryMap == "TRUE" ]; then
# 				ofile="$odir/cortmap_${metric}_resamplesummary_AB${yearsBeforeAB}nTime${nTime}_${analysname}.png"
# 			else
# 				ofile="$odir/${metric}_cortmapAB${yearsBeforeAB}_nTime${nTime}_${analysname}.png"
# 			fi
# 			# if [ ! -e "$ofile" ]; then
# 				if [ $onClust == 1 ]; then
# 					cmd="sbatch X-makeCortmap.sh"
# 				else
# 					cmd="sh X-makeCortmap.sh"
# 				fi
# 				$cmd $analysname $yearsBeforeAB $metric $nTime $summaryMap $empiricalMap
# 			# else
# 			# 	echo "$ofile exists. Skipping."
# 			# fi
# 		done
# 	done
# fi


# if [ $step = 6 ]; then
# # #other maps
# 	echo "----------- making other ggseg maps (yearsBeforeAB) -----------"
# 	summaryMap="FALSE"
# 	empiricalMap="FALSE"
# 	omaps="Yintonly Yintonlycentcor noage noagecentcor Ymanchangenoout Ymanchangenooutcentcor \
# 	Yslopoys Yslopoyscentcor Yslointooys Yslointooyscentcor"

# 	logd="logs-06"
# 	if [ ! -d $logd ]; then 
# 		mkdir $logd
# 	fi
# 	for yearsBeforeAB in $yearcuts; do
# 		for omap in $omaps; do
# 			ifile="${metric}_cortmapAB${yearsBeforeAB}_nTime${nTime}_${analysname}_${omap}"
# 			ofile="$odir/${ifile}.png"
# 			mapName="${ifile}.csv"
# 			if [ ! -e "$ofile" ]; then
# 				if [ $onClust == 1 ]; then
# 					cmd="sbatch XX-makeCortmap.sh"
# 				else
# 					cmd="sh XX-makeCortmap.sh"
# 				fi
# 			else
# 				echo "$ofile exists. Skipping."
# 			fi
# 			$cmd $analysname $yearsBeforeAB $metric $nTime $summaryMap $empiricalMap $mapName
# 		done
# 	done
# fi


# if [ $step = 7 ]; then
# 	echo "----------- making combined figures (yearsBeforeAB) -----------"

# 	odir="$resdir/mapPlots/$analysname/$metric"
# 	fdir="$base/figures/$analysname/$metric"
# 	# echo $odir; exit 0
# 	echo $odir
# 	echo $metric

# 	logd="logs-07"
# 	if [ ! -d $logd ]; then 
# 		mkdir $logd
# 	fi
# 	# if [ ! -d "$fdir" ]; then

# 		if [ $onClust == 1 ]; then
# 			cmd="sbatch X-create_figures.sh"
# 		else
# 			cmd="sh X-create_figures.sh"
# 		fi
# 		$cmd $odir $metric
# 	# else
# 		# echo "$fdir exists. Skipping."
# 	# fi

# fi


# if [ $step = 8 ]; then
# 	echo "----------- making combined figures (yearsBeforeAB) -----------"

# 	odir="$resdir/mapPlots/$analysname/$metric"
# 	fdir="$base/figures/$analysname/$metric"
# 	# echo "OUTPUT_DIR='${odir}'"; echo 'omap="Yintonly"'; echo 'metric="thickness"'; exit 0
# 	echo $odir
# 	echo $metric
# 	omaps="Yintonly noage Ymanchangenoout Yslopoys Yslointooys" #note that centcor versions get detected
	
# 	logd="logs-08"
# 	if [ ! -d $logd ]; then 
# 		mkdir $logd
# 	fi
# 	for omap in $omaps; do
			
# 		# if [ ! -d "$fdir" ]; then

# 			if [ $onClust == 1 ]; then
# 				cmd="sbatch XX-create_figures.sh"
# 			else
# 				cmd="sh XX-create_figures.sh"
# 			fi
# 			$cmd $odir $metric $omap
# 		# else
# 			# echo "$fdir exists. Skipping."
# 		# fi
# 	done
# fi


if [ $step = 9 ]; then
	exportmenow() {
		cp -R $@ /tsd/p274/data/durable/file-export
	}
	figdir=$base/figures
	
	ls $figdir/$analysname/$metric
	ls $figdir/$analysname/$metric | wc -l; sleep 3
	cd $figdir
	echo "> zipping figures ${analysname}.tgz"; sleep 1
	tar czvf ${analysname}.tgz ${analysname}
	echo "> exporting figures ${analysname}.tgz"
	exportmenow ${analysname}.tgz
	cd $scrdir
fi


if [ $step = 10 ]; then
	echo "----------- running pPerm -----------"

	regions="lh_bankssts_thickness.aparcnative71 \
	lh_caudalanteriorcingulate_thickness.aparcnative71 \
	lh_caudalmiddlefrontal_thickness.aparcnative71 \
	lh_cuneus_thickness.aparcnative71 \
	lh_entorhinal_thickness.aparcnative71 \
	lh_fusiform_thickness.aparcnative71 \
	lh_inferiorparietal_thickness.aparcnative71 \
	lh_inferiortemporal_thickness.aparcnative71 \
	lh_isthmuscingulate_thickness.aparcnative71 \
	lh_lateraloccipital_thickness.aparcnative71 \
	lh_lateralorbitofrontal_thickness.aparcnative71 \
	lh_lingual_thickness.aparcnative71 \
	lh_medialorbitofrontal_thickness.aparcnative71 \
	lh_middletemporal_thickness.aparcnative71 \
	lh_parahippocampal_thickness.aparcnative71 \
	lh_paracentral_thickness.aparcnative71 \
	lh_parsopercularis_thickness.aparcnative71 \
	lh_parsorbitalis_thickness.aparcnative71 \
	lh_parstriangularis_thickness.aparcnative71 \
	lh_pericalcarine_thickness.aparcnative71 \
	lh_postcentral_thickness.aparcnative71 \
	lh_posteriorcingulate_thickness.aparcnative71 \
	lh_precentral_thickness.aparcnative71 \
	lh_precuneus_thickness.aparcnative71 \
	lh_rostralanteriorcingulate_thickness.aparcnative71 \
	lh_rostralmiddlefrontal_thickness.aparcnative71 \
	lh_superiorfrontal_thickness.aparcnative71 \
	lh_superiorparietal_thickness.aparcnative71 \
	lh_superiortemporal_thickness.aparcnative71 \
	lh_supramarginal_thickness.aparcnative71 \
	lh_frontalpole_thickness.aparcnative71 \
	lh_temporalpole_thickness.aparcnative71 \
	lh_transversetemporal_thickness.aparcnative71 \
	lh_insula_thickness.aparcnative71 \
	lh_MeanThickness_thickness.aparcnative71 \
	rh_bankssts_thickness.aparcnative71 \
	rh_caudalanteriorcingulate_thickness.aparcnative71 \
	rh_caudalmiddlefrontal_thickness.aparcnative71 \
	rh_cuneus_thickness.aparcnative71 \
	rh_entorhinal_thickness.aparcnative71 \
	rh_fusiform_thickness.aparcnative71 \
	rh_inferiorparietal_thickness.aparcnative71 \
	rh_inferiortemporal_thickness.aparcnative71 \
	rh_isthmuscingulate_thickness.aparcnative71 \
	rh_lateraloccipital_thickness.aparcnative71 \
	rh_lateralorbitofrontal_thickness.aparcnative71 \
	rh_lingual_thickness.aparcnative71 \
	rh_medialorbitofrontal_thickness.aparcnative71 \
	rh_middletemporal_thickness.aparcnative71 \
	rh_parahippocampal_thickness.aparcnative71 \
	rh_paracentral_thickness.aparcnative71 \
	rh_parsopercularis_thickness.aparcnative71 \
	rh_parsorbitalis_thickness.aparcnative71 \
	rh_parstriangularis_thickness.aparcnative71 \
	rh_pericalcarine_thickness.aparcnative71 \
	rh_postcentral_thickness.aparcnative71 \
	rh_posteriorcingulate_thickness.aparcnative71 \
	rh_precentral_thickness.aparcnative71 \
	rh_precuneus_thickness.aparcnative71 \
	rh_rostralanteriorcingulate_thickness.aparcnative71 \
	rh_rostralmiddlefrontal_thickness.aparcnative71 \
	rh_superiorfrontal_thickness.aparcnative71 \
	rh_superiorparietal_thickness.aparcnative71 \
	rh_superiortemporal_thickness.aparcnative71 \
	rh_supramarginal_thickness.aparcnative71 \
	rh_frontalpole_thickness.aparcnative71 \
	rh_temporalpole_thickness.aparcnative71 \
	rh_transversetemporal_thickness.aparcnative71 \
	rh_insula_thickness.aparcnative71 \
	rh_MeanThickness_thickness.aparcnative71"
	# pcaL \
	# pcaR"
	regionsarr=($regions)

	for yearsBeforeAB in $(seq 1 10); do

		for roi in $regions; do
			odir="$resdir/cortmaps_resample/$analysname/${metric}_nullmod/yearsBeforeAB${yearsBeforeAB}nTime${nTime}/${roi}"
			ofile1="$odir/nullresult_yearsBeforeAB${yearsBeforeAB}nTime${nTime}-${roi}.csv"
			ofile2="$odir/nullresultcentcor_yearsBeforeAB${yearsBeforeAB}nTime${nTime}-${roi}.csv"
			#echo $ofile1; echo $ofile2
			
			ofile3="$odir/region_yearsBeforeAB${yearsBeforeAB}nTime${nTime}-${roi}.csv"
			
			
			if [[ ! -e $ofile1 ]] || [[ ! -e $ofile2 ]]; then
				if [[ -e $ofile3 ]]; then
					echo -n "Submitting job #$i... "
					if [ $onClust == 1 ]; then
						cmd="sbatch 06-regionalWildBootstrap.batch"
					else
						cmd="sh 06-regionalWildBootstrap.batch"
					fi
					$cmd $roi $yearsBeforeAB $analysname $metric $nTime
				else
					echo "$ofile3 does not exist. Will not run."
				fi
			else
				echo "at least one outfile exists: $ofile1. Skipping."
			fi
		done
	done
fi


if [ $step = 11 ]; then
	echo "----------- creating empirical p-value map -----------"
	for yearsBeforeAB in $(seq 1 10); do
		Rscript 07-empiricalPmaps.r $yearsBeforeAB $analysname $metric $nTime
	done
fi


# if [ $step = 12 ]; then
# 	echo "----------- making ggseg maps (empirical yearsBeforeAB) -----------"
	
# 	odir="$resdir/mapPlots/$analysname/${metric}_null"
# 	empiricalMap=TRUE
# 	for yearsBeforeAB in $(seq 1 10); do
# 		for summaryMap in FALSE; do
# 			ofile="$odir/empirical_${metric}_cortmap_AB${yearsBeforeAB}_nTime${nTime}_${analysname}.png"
			
# 			if [ ! -e "$ofile" ]; then
# 				if [ $onClust == 1 ]; then
# 					cmd="sbatch X-makeCortmap.sh"
# 				else
# 					cmd="sh X-makeCortmap.sh"
# 				fi
# 				$cmd $analysname $yearsBeforeAB $metric $nTime $summaryMap $empiricalMap
# 			else
# 				echo "$ofile exists. Skipping."
# 			fi
# 		done
# 	done
# fi

