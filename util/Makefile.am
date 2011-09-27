

if FULLDIST
bin_SCRIPTS = gmap_compress gmap_uncompress gmap_process gmap_setup gmap_build gmap_reassemble md_coords fa_coords psl_splicesites psl_splices dbsnp_iit sam_merge sam_restore
else
bin_SCRIPTS = gmap_compress gmap_uncompress gmap_process gmap_setup gmap_build gmap_reassemble md_coords fa_coords psl_splicesites psl_splices dbsnp_iit sam_merge sam_restore
endif

gmap_compress: gmap_compress.pl
	cp gmap_compress.pl gmap_compress
	chmod +x gmap_compress

gmap_uncompress: gmap_uncompress.pl
	cp gmap_uncompress.pl gmap_uncompress
	chmod +x gmap_uncompress

gmap_process: gmap_process.pl
	cp gmap_process.pl gmap_process
	chmod +x gmap_process

gmap_setup: gmap_setup.pl
	cp gmap_setup.pl gmap_setup
	chmod +x gmap_setup

gmap_build: gmap_build.pl
	cp gmap_build.pl gmap_build
	chmod +x gmap_build

gmap_update: gmap_update.pl
	cp gmap_update.pl gmap_update
	chmod +x gmap_update

gmap_reassemble: gmap_reassemble.pl
	cp gmap_reassemble.pl gmap_reassemble
	chmod +x gmap_reassemble

md_coords: md_coords.pl
	cp md_coords.pl md_coords
	chmod +x md_coords

fa_coords: fa_coords.pl
	cp fa_coords.pl fa_coords
	chmod +x fa_coords

psl_splicesites: psl_splicesites.pl
	cp psl_splicesites.pl psl_splicesites
	chmod +x psl_splicesites

psl_splices: psl_splices.pl
	cp psl_splices.pl psl_splices
	chmod +x psl_splices

dbsnp_iit: dbsnp_iit.pl
	cp dbsnp_iit.pl dbsnp_iit
	chmod +x dbsnp_iit

sam_merge: sam_merge.pl
	cp sam_merge.pl sam_merge
	chmod +x sam_merge

sam_restore: sam_restore.pl
	cp sam_restore.pl sam_restore
	chmod +x sam_restore


if FULLDIST
CLEANFILES = gmap_compress gmap_uncompress gmap_process gmap_setup gmap_build gmap_reassemble md_coords fa_coords psl_splicesites psl_splices dbsnp_iit sam_merge sam_restore
else
CLEANFILES = gmap_compress gmap_uncompress gmap_process gmap_setup gmap_build gmap_reassemble md_coords fa_coords psl_splicesites psl_splices dbsnp_iit sam_merge sam_restore
endif

