

if FULLDIST
bin_SCRIPTS = gmap_compress gmap_uncompress gmap_process gmap_setup gmap_update md_coords fa_coords
else
bin_SCRIPTS = gmap_compress gmap_uncompress gmap_process gmap_setup gmap_update md_coords fa_coords
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

gmap_update: gmap_update.pl
	cp gmap_update.pl gmap_update
	chmod +x gmap_update

md_coords: md_coords.pl
	cp md_coords.pl md_coords
	chmod +x md_coords

fa_coords: fa_coords.pl
	cp fa_coords.pl fa_coords
	chmod +x fa_coords

if FULLDIST
CLEANFILES = gmap_compress gmap_uncompress gmap_process gmap_setup gmap_update md_coords fa_coords
else
CLEANFILES = gmap_compress gmap_uncompress gmap_process gmap_setup gmap_update md_coords fa_coords
endif

