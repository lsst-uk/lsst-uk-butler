#! /bin/bash
# setup skymap
python src/autoIngest.py -s

# setup the butler
python src/autoIngest.py -b

# ingest raw data
python src/autoIngest.py -w

# register and ingest coadds
python src/autoIngest.py -r c
python src/autoIngest.py -i c

# define visits
python src/autoIngest.py -v

# register and ingest image data
python src/autoIngest.py -r i
python src/autoIngest.py -i i

# ingest calib data
python src/autoIngest.py -c

# ingest reference monster catalogues
python src/autoIngest.py -m

# register and ingest all other data
python src/autoIngest.py -r a
python src/autoIngest.py -i a
