# load useful packages ----
library(tidyverse)
round0 = function(x,z){
	round(x/z,0)*z
}


# Convert EDFs ----

#define a function to convert a given edf
convert_edfs = function(file,.pb=NULL){
	if ((!is.null(.pb)) && inherits(.pb, "Progress") && (.pb$i < .pb$n)) .pb$tick()$print()
	system(
		command = paste("./edf2asc",file)
	)
	system(
		command = paste("rm",file)
	)
	file = gsub('edf','asc',file)
	system(
		command = paste("gzip --fast",file)
	)
}

#get list of files to convert
files = list.files(
	path = '../_Data'
	, pattern = '.edf'
	, full.names = T
	, recursive = T
)

#run walk through files, with progress bar
pb = dplyr::progress_estimated(length(files))
purrr::walk(
	.x = files
	, .f = convert_edfs
	, .pb = pb
)

# Extract info from ASCs ----

#define a function to process a given asc
pre_process_eye_data = function(x,.pb=NULL){
	if ((!is.null(.pb)) && inherits(.pb, "Progress") && (.pb$i < .pb$n)) .pb$tick()$print()
	if(!file.exists(gsub('.asc.gz','_blinks.txt.gz',x))){
		#read everything
		a = readLines(x)

		#remove header
		a = a[substr(a,1,2)!="**"]

		#remove calibration headers
		a = a[substr(a,1,7)!=">>>>>>>"]

		#remove blank lines
		a = a[a!=""]

		#remove lines that start with blanks (from clibration)
		a = a[substr(a,1,1)!=" "]
		a = a[substr(a,1,1)!="\t"]

		#isolate & write msgs
		subset = substr(a,1,3)=="MSG"
		msgs = a[subset]
		write(msgs,file = gsub('.asc.gz','_msgs.txt',x))
		system(
			command = paste("gzip --fast",gsub('.asc.gz','_msgs.txt',x))
		)
		rm(msgs)
		a = a[!subset]
		gc()

		#isolate & write inputs
		subset = substr(a,1,5)=="INPUT"
		inputs = a[subset]
		write(inputs,file = gsub('.asc.gz','_inputs.txt',x))
		system(
			command = paste("gzip --fast",gsub('.asc.gz','_inputs.txt',x))
		)
		rm(inputs)
		a = a[!subset]
		gc()

		#isolate & write samples
		subset = substr(a,nchar(a)-2,nchar(a))=="..."
		samples = a[subset]
		write(samples,file = gsub('.asc.gz','_samples.txt',x))
		system(
			command = paste("gzip --fast",gsub('.asc.gz','_samples.txt',x))
		)
		rm(samples)
		a = a[!subset]
		gc()

		#isolate & write saccades
		subset = substr(a,1,5)=="ESACC"
		saccades = a[subset]
		write(saccades,file = gsub('.asc.gz','_saccades.txt',x))
		system(
			command = paste("gzip --fast",gsub('.asc.gz','_saccades.txt',x))
		)
		rm(saccades)
		a = a[!subset]
		gc()

		#isolate & write blinks
		subset = substr(a,1,6)=="EBLINK"
		temp = a[subset]
		write(temp,file = gsub('.asc.gz','_blinks.txt',x))
		system(
			command = paste("gzip --fast",gsub('.asc.gz','_blinks.txt',x))
		)
		rm(temp)
		a = a[!subset]
		gc()
	}

}

#get list of files to process
files = list.files(
	path = '../_Data'
	, pattern = '.asc.gz'
	, full.names = T
	, recursive = T
)

#walk through files, with progress bar
pb = dplyr::progress_estimated(length(files))
purrr::walk(
	.x = files
	, .f = pre_process_eye_data
	, .pb = pb
)

# Get data ----

#define a function to process a given subject's data
get_subject_data = function(folder,.pb=NULL){
	if ((!is.null(.pb)) && inherits(.pb, "Progress") && (.pb$i < .pb$n)) .pb$tick()$print()

	#get samples
	sample_file = list.files(folder,full.names=T,recursive = T,pattern='samples')
	samples  = read.table(sample_file,stringsAsFactors = F,na.strings = '.')[,1:4]
	names(samples) = c('time','x','y','pupil')
	samples$x = samples$x-1024/2
	samples$y = samples$y-768/2
	samples$dist = sqrt( (samples$x^2) + (samples$y^2) )
	samples$out_of_box = (abs(samples$x)>100) | (abs(samples$y)>100)
	# mean(samples$out_of_box,na.rm=T)
	# mean(samples$dist>100,na.rm=T)

	#get blinks
	blink_file = list.files(folder,full.names=T,recursive = T,pattern='blinks')
	blink_times = read.table(blink_file)[,3]

	# #get saccades
	# saccade_file = list.files(file,full.names=T,recursive = T,pattern='saccades')
	# saccades = read.table(saccade_file)[,c(3,8,9)]
	# names(saccades) = c('start','x','y')

	#get msgs
	msg_file = list.files(folder,full.names=T,recursive = T,pattern='msgs')
	msgs  = readLines(msg_file)
	subset = grepl('StudyWord',msgs)
	subset = which(subset)

	words = msgs[subset+1]
	words = gsub('MSG\t','',words)
	words = tibble::as_tibble(words)
	words %>%
		tidyr::separate(
			value
			, into = c('word_time','word')
			, convert = T
		) -> words


	instructions = msgs[subset+5]
	instructions = gsub('MSG\t','',instructions)
	instructions = tibble::as_tibble(instructions)
	instructions %>%
		tidyr::separate(
			value
			, into = c('instruction_time','instruction')
			, convert = T
		) -> instructions

	msgs = bind_cols(words,instructions)
	rm(words)
	rm(instructions)

	#check that the diff is ~2s
	# summary(msgs$instruction_time - msgs$word_time)

	out = list()
	for(i in 1:nrow(msgs)){
		t0 = msgs$instruction_time[i]
		critical_blink = any( ((blink_times-t0)<2000) & ((blink_times-t0)>(-500)) )
		#	if(!critical_blink){ #no blinks during this period
		subset = ((samples$time-t0)<2000) & ((samples$time-t0)>(-500))
		samps = samples[subset,]
		samps = samps[!is.na(samps$dist),]
		# samps = samps[samps$dist<100,]
		if(nrow(samps)>0){
			samps$time = samps$time-t0
			#samps$time = round0(samps$time,1)
			samps %>%
				group_by(time) %>%
				summarise(
					pupil = mean(pupil)
					, x = mean(x)
					, y = mean(y)
					, dist = mean(dist)
				) -> samps
			samps$pupil = samps$pupil-mean(samps$pupil[samps$time<0])
			out[[length(out)+1]] = tibble::as_tibble(data.frame(
				trial = i
				, critical_blink = critical_blink
				, word = msgs$word[i]
				, instruction = msgs$instruction[i]
				, time = samps$time
				, pupil = samps$pupil
				, x = samps$x
				, y = samps$y
				, dist = samps$dist
			))
		}
		#	}
	}
	out2 = dplyr::bind_rows(out)
	rm(out)
	gc()
	out2 %>%
		dplyr::mutate(
			stim = case_when(
				.$instruction %in% c('cr','cf') ~ 'non-word'
				, .$instruction %in% c('r','f') ~ 'word'
			)
			, instruction = case_when(
				.$instruction %in% c('cr','r') ~ 'remember'
				, .$instruction %in% c('cf','f') ~ 'forget'
			)
		) -> out2
	out2$id = strsplit(folder,'/')[[1]][4]
	return(out2)
}

#get list of folders to process
folders = list.files(
	path = '../_Data'
	, full.names = T
)

#remember to exclude dfp106 (glasses) and dfp107 (crash)
folders = folders[!grepl("102tnt",tolower(folders))]
folders = folders[!grepl("105tnt",tolower(folders))]


#collect data from each folder, with progress bar
pb = dplyr::progress_estimated(length(folders))
a = purrr::map_df(
	.x = folders
	, .f = get_subject_data
	, .pb = pb
)

hist(a$dist,br=100)
mean(a$dist>100)
save(a,file='a.rdata')
#load('a.rdata')

a %>%
	dplyr::group_by(
		instruction
		, stim
		, id
	) %>%
	dplyr::summarise(
		value = mean(critical_blink,na.rm=T)
		, count = n()
	) %>%
	ggplot(
		mapping = aes(
			x = instruction
			,  y = value
		)
	)+
	geom_boxplot()+
	geom_point(alpha=.5)+
	facet_grid(.~stim)+
	theme(
		aspect.ratio = 1
	)

a %>%
	ggplot()+
	geom_smooth(
		# geom_line(
		mapping = aes(
			x = time
			, y = dist
			, colour = instruction
		)
	)+
	facet_grid(.~stim)+
	theme(
		aspect.ratio = 1
	)

