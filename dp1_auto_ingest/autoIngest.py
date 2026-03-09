import getopt
import os
import sys
import glob
import copy
import subprocess
import re
import traceback
import pickle
import io
import selectors
import shlex

from astropy.table import Table
import pandas as pd
from collections import defaultdict, namedtuple, OrderedDict
from dotenv import load_dotenv
from contextlib import redirect_stdout
from copy import deepcopy
#from subprocess import Popen, PIPE, STDOUT

from lsst.daf.butler import Butler
#import lsst.daf.butler as Butler

#-----------------------------------------------------------------------------
dotenvfile = '.env'
#dotenvfile = '.dp1testenv'
load_dotenv(os.path.join(os.environ['HOME'], dotenvfile))
#
# data directories
#
repo = os.environ['REPO']
baseDir = os.environ['DATA']
defDataDir = 'dp1/LSSTComCam'
defRawDataDir =  'raw/LSSTComCam'
defRefDir = 'raw/refcats'
defSkyMapDir = 'dp1/skymaps/skyMap'
defCalibDir = 'LSSTComCam/calib'
defRawDir = 'LSSTComCam/raw'
csvBaseDir = os.path.join(repo, 'auto_csv_files')
ecsvBaseDir = os.path.join(repo, 'auto_ecsv')
defSrcDir = os.path.join(repo, 'src')
butlerPath = '/data/sandbox/lsst_stack_29/conda/envs/lsst-scipipe-10.1.0/share/eups/Linux64/daf_butler/g6dd59efbe6+ec144cf465/bin/butler'
#
# general lsst variables
#
skymap = 'lsst_cells_v1'
instrument = 'LSSTComCam'
release = 'DP1'
defCollection = f'{instrument}/{release}'
fullinstrument = 'lsst.obs.lsst.LsstComCam'
detectorIdComCam = {'R22_S00': 0, 'R22_S01': 1, 'R22_S02': 2,
                    'R22_S10': 3, 'R22_S11': 4, 'R22_S12': 5,
                    'R22_S20': 6, 'R22_S21': 7, 'R22_S22': 8}
firstIngDsTypes = [['deep_coadd', 'template_coadd'],
                   ['visit_image', 'difference_image']]
#
# data files
#
dataTypesFile = os.path.join(defSrcDir, 'allDataTypesUSDF.list')
skyMapPickleFile = os.path.join(
    baseDir, defSkyMapDir, f'skyMap_{skymap}_skymaps.pickle')
skyMapConfigFile = os.path.join(defSrcDir, f'skyMap_{skymap}_skymaps.config')
refCats = ['the_monster_20250219']
#
# support variables
#
list_of_dims = ['new_url', 'instrument', 'skymap', 'band', 'patch',
                'tract', 'visit', 'physical_filter', 'day_obs', 'detector',
                'htm7', 'exposure']
dim_defs = [''] * len(list_of_dims)
DefaultDimension = OrderedDict(zip(list_of_dims, dim_defs))
verbose = 0
#
# set panda options for table views
#
pd.options.display.max_rows = None
pd.options.display.max_columns = None
pd.options.display.max_colwidth = None
pd.options.display.expand_frame_repr = False
#=============================================================================
#-----------------------------------------------------------------------------

def theButler(repository, collection=None):
    if verbose > 5:
        print(">>>",repository)
    butler = Butler(repository, collections=collection)
    assert butler is not None
    return butler

#-----------------------------------------------------------------------------

def butler_setup(testRun=False):

    # setup the butler
    butlyaml = os.path.join(repo, 'butler.yaml')
    cmd = f'{butlerPath} create --seed-config {repo}/butler-seed.yaml --override {repo}'
    cmdargs =  shlex.split(cmd)
    if testRun:
        print(f'Deleting and recreating {butlyaml}')
        print(f'>> {cmd}')
    else:
        # delete old file
        if os.path.exists(butlyaml):
            os.remove(butlyaml)
            print(f'Deleted {butlyaml}')

        # recreate butler.yaml
        result = run_shell_cmd(cmdargs)
        print(f'Recreated {butlyaml}')
    
    # register instrument
    cmd = f"{butlerPath} register-instrument {repo} '{fullinstrument}'"
    cmdargs =  shlex.split(cmd)
    if testRun:
        print(f'>> {cmd}')
    else:
        print(f'>> {cmdargs}')
        result = run_shell_cmd(cmdargs)
        print(f'Registered {fullinstrument}')

    # register skymap
    if os.path.exists(skyMapConfigFile):
        cmd = f'{butlerPath} register-skymap -C {skyMapConfigFile} {repo}'
        cmdargs =  shlex.split(cmd)
        if testRun:
            print(f'>> {cmd}')
        else:
            print('Registering skymap....')
            print(f'>> {cmdargs}')
            result = run_shell_cmd(cmdargs)
            print(f'Skymap {skyMapConfigFile} registered')
    else:
        print('[ERROR] Skymap config file {skyMapConfigFile} not found!')
        raise SystemExit

#-----------------------------------------------------------------------------

def butler_ingest_raw(testRun=False):

    # ingest raw data
    rawdir = os.path.join(baseDir, defRawDataDir)
    ingestmethod = '--transfer symlink --processes 4'
    if os.path.exists(rawdir):
        cmd = f'{butlerPath} ingest-raws {ingestmethod} {repo} {rawdir}'
        cmdargs =  shlex.split(cmd)
        if testRun:
            print(f'>> {cmd}')
        else:
            print('Ingesting raws....')
            print(f'>> {cmdargs}')
            result = run_shell_cmd(cmdargs)
            print("Raws ingested")
    else:
        print(f'[ERROR] Raw path {rawdir} not found!') 

#-----------------------------------------------------------------------------

def butler_define_visits(testRun=False):
    # define visits
    collection = os.path.join(defRawDir, 'all')
    cmd = f'butler define-visits {repo} {instrument} --collections {collection}'

    cmdargs =  shlex.split(cmd)
    if testRun:
        print(f'>> {cmd}')
    else:
        print("Defining visits...")
        print(f'>> {cmdargs}')
        result = run_shell_cmd(cmdargs)
        print("Visits defined")
        
#-----------------------------------------------------------------------------

def butler_query_one_dataset(butler, datatype, testRun=False):
    cmd = f"butler.registry.queryDatasetTypes('{datatype}')"
    print(f'>> {cmd}')
    thedatasettype = ''
    if not testRun:
        thedatasettype = butler.registry.queryDatasetTypes(datatype)
        print(f'>> {thedatasettype}')
    return True if thedatasettype else False

#-----------------------------------------------------------------------------

def butler_register_one_dataset(datatype, theDataSetType, testRun=False):
    dimensions, format = theDataSetType
    print(f'Registering {datatype}')
    dims = ' '.join(dimensions)

    cmd = f'butler register-dataset-type {repo} {datatype} {format} {dims}'
    cmdargs =  shlex.split(cmd)
    if testRun:
        print(f'>> {cmd}')
        result = 0
    else:
        print(f'Registering {datatype}')
        print(f'>> {cmdargs}')
        result = run_shell_cmd(cmdargs)
    return result

#-----------------------------------------------------------------------------

def butler_ingest_one_dataset(datatype, testRun=False):
    print(f'Ingesting {datatype}')
    ecsvfile =  os.path.join(ecsvBaseDir, f'{datatype}.ecsv')

    # get the 'runpath'
    with open(ecsvfile, 'r') as f:
        thelines = f.readlines()
    thefiles = [line for line in thelines if line.startswith(baseDir)]
    dp1Dir = os.path.join(baseDir, 'dp1/')
    comdirs = os.path.commonpath(thefiles).replace(dp1Dir,'').split('/')
    theIdx = comdirs.index(next(filter(lambda dir: dir.startswith('DM-'),
                                       comdirs), None))
    if comdirs[theIdx-1].startswith('v29_'):
        maxIdx = theIdx+1
    else:
        maxIdx = theIdx
    runpath = '/'.join(comdirs[:maxIdx+1])

    cmd = f'butler ingest-files {repo} {datatype} {runpath} {ecsvfile}'
    cmdargs =  shlex.split(cmd)
    if testRun:
        print(f'>> {cmd}')
        result = 0
    else:
        print(f'Ingesting {datatype}')
        print(f'>> {cmdargs}')
        result = run_shell_cmd(cmdargs)
    return result

#-----------------------------------------------------------------------------

def butler_reg_ing_datasets(dataSetTypes, testRun=False):
    for datatype in dataSetTypes:
        print(f'Registering & Ingesting {datatype}:')
        # check if registered
        butler = theButler(repo)
        exists = butler_query_one_dataset(butler, datatype, False)
        if exists:
            print(f'{datatype} already exists!')
            testRun = True

        # register
        res = butler_register_one_dataset(
            datatype, dataSetTypes[datatype], testRun)

        # ingest
        res = butler_ingest_one_dataset(datatype, testRun)

#-----------------------------------------------------------------------------

def butler_register_datasets(dataSetTypes, testRun=False):
    print(f'Registering {dataSetTypes.keys()}')
    for datatype in dataSetTypes:
        butler_register_one_dataset(datatype, theDataSetTypes[datatype],
                                    testRun)
    print(f'Registered all datasetTypes {dataSetTypes.keys()}')

#-----------------------------------------------------------------------------

def butler_ingest_datasets(dataSetTypes, testRun=False):
    print(f'Ingesting {dataSetTypes.keys()}')
    for datatype in dataSetTypes:
        butler_ingest_one_dataset(datatype, testRun)
    print(f'Ingested all datasetTypes {dataSetTypes.keys()}')

#-----------------------------------------------------------------------------

def butler_ingest_calib(testRun=False):
    calibdirs = os.listdir(os.path.join(baseDir, defDataDir, 'calib'))
    print(f'Ingesting calibrations...')
    for calib in calibdirs:
        # ingest calib data
        if calib.startswith('DM'):
            cmd = f'butler write-curated-calibrations {repo} {fullinstrument} --label {calib}'
            cmdargs =  shlex.split(cmd)
            if testRun:
                print(f'>> {cmd}')
            else:
                print(f'>> {cmdargs}')
                result = run_shell_cmd(cmdargs)

            print('Creating collection chain...')
            # create collection chain
            coll1 = [os.path.join(defCalibDir, x)
                     for x in calibdirs if x.startswith('DM')]
            coll2 = [os.path.join(x, 'unbounded') for x in coll1]
            collections = ','.join(coll1 + coll2)
            cmd = f'butler collection-chain {repo} {defCalibDir} {collections}'
            cmdargs =  shlex.split(cmd)
            if testRun:
                print(f'>> {cmd}')
            else:
                print(f'>> {cmdargs}')
                result = run_shell_cmd(cmdargs)
    print('Calibrations ingested and chains created')
            
#-----------------------------------------------------------------------------

def butler_ingest_refcats(testRun=False):
    refdir = os.path.join(baseDir, defRefDir)
    # ingest reference catalogues
    for theRefcat in refCats:
        ecsvfile =  os.path.join(ecsvBaseDir, f'{theRefcat}.ecsv')
        refCatsDir = os.path.basename(defRefDir)
        runpath = os.path.join(refCatsDir, theRefcat)
        prefix = os.path.join(refdir, theRefcat)
        cmd = f'butler ingest-files {repo} {theRefcat} {runpath} --prefix {prefix}/ --transfer symlink  {ecsvfile}'
        cmdargs =  shlex.split(cmd)
        if testRun:
            print(f'>> {cmd}')
        else:
            print(f'Ingesting {theRefcat}...')
            print(f'>> {cmdargs}')
            result = run_shell_cmd(cmdargs)

        # create collection chain
        print(f'Creating {theRefcat} collection chain...')
        cmd = f'butler collection-chain {repo} {refCatsDir} {runpath}'
        cmdargs =  shlex.split(cmd)
        if testRun:
            print(f'>> {cmd}')
        else:
            print(f'>> {cmdargs}')
            result = run_shell_cmd(cmdargs)
        print(f'{theRefcat} ingested and chains created')

#-----------------------------------------------------------------------------

def butler_query(mode=['i'], dataSetTypes=[], testRun=False):
    # query what's ingested
    butler = theButler(repo)#, defCollection)
    if 'i' in mode:
        # information about main release collection
        print(f'The */* collection:')
        cmd = "butler.collections.query_info('*/*')"
        print(f'>> {cmd}')
        if not testRun:
            thecollection = butler.collections.query_info('*/*')
            for coll in thecollection:
                print(f'   {coll}')

    if 's' in mode:
        # get skymap info
        print('Skymap info:')
        cmd = "butler.query_dimension_records('skymap')"
        print(f'>> {cmd}')
        if not testRun:
            dataset_type = butler.get_dataset_type('skyMap')
            for dimension in dataset_type.dimensions.data_coordinate_keys:
                print('dimension = ', dimension)
                print(butler.dimensions[dimension].schema)
                print(' ')
            
            info = butler.query_dimension_records('skymap')
            print(info)

    if 'c' in mode:
        # information about all collections
        print('All collections:')
        cmd = "butler.collections.query('*')"
        print(f'>> {cmd}')
        if not testRun:
            all_collections = butler.collections.query('*')
            for collection in all_collections:
                print(f'>> {collection}')
                
    if 't' in mode:
        # information about all dataset types
        print('All dataset types:')
        cmd = "butler.registry.queryDatasetTypes()"
        print(f'>> {cmd}')
        if not testRun:
            all_datasettypes = butler.registry.queryDatasetTypes()
            for dst in all_datasettypes:
                print(f'>> {dst}')

    if 'd' in mode:
        # information if dataset contains data
        print('Data ingested check:')
        cmd = "butler.query_datasets(dst, limit=10)"
        print(f'>> {cmd}')
        if not testRun:
            all_datasettypes = butler.registry.queryDatasetTypes()
            for dst in all_datasettypes:
                dataset_refs = butler.query_datasets(dst, limit=10)
                print(f'>> {dst}: {len(dataset_refs)}')
                del dataset_refs

#-----------------------------------------------------------------------------

def change_base_directory(url, dstype, dstype_dims):
    directories=url.split('/')[5:]
    filename = url.rsplit('/',1)[1]
    new_url = os.path.join(baseDir, *directories)

    if 'raw' in directories[0]:
        dtype = directories[2]
    elif 'DM' in directories[5] and 'runs' in directories[2]:
        dtype = directories[6]
    elif 'v29' in directories[5]:
        dtype = directories[8]
    elif 'DM' in directories[4]:
        if 'fgcm' in directories[5]:
            dtype = directories[5]
        elif 'standard' in directories[5]:
            dtype = directories[6]
    elif 'DM' in directories[3]:
        if any(x in directories[5]
               for x in ['flat', 'bfk', 'defect', 'cti', 'ptc', 'dark',
                         'linearizer', 'bias']):
            dtype = directories[7]
        elif 'unbounded' in directories[4]:
            dtype = directories[5]
        elif 'curated' in directories[4]:
            dtype = directories[6]

    if verbose > 2:
        print("????",dstype,"::",directories)
        print("!!!!!",dtype,"::",directories,"####",dstype)

    if dtype != dstype:
        print(f'<ERROR> DataSetTypes differ: given {dstype} <> read {dtype}')
        raise SystemExit
    
    if verbose > 2:
        print(directories,':::',dtype)

    Dimension = deepcopy(DefaultDimension)
    Dimension.update([('new_url', new_url),
                      ('instrument', instrument),
                      ('skymap', skymap)])
    dst_dims = [dim for dim in dstype_dims
                if dim not in ['instrument', 'skymap']]
    if verbose > 3:
        print(Dimension)

    filenamedata = filename.replace(dtype, '')
    if verbose > 2:
        print("<##>", filenamedata)
        print("<><>", dst_dims)
        
    if set(dst_dims) == set(['band', 'day_obs', 'physical_filter',
                             'tract', 'patch', 'visit']):
        Dimension.update([
            ('visit', directories[-2]),
            ('physical_filter', directories[-3]),
            ('band', directories[-4]),
            ('day_obs', directories[-5]),
            ('patch', directories[-6]),
            ('tract', directories[-7])
        ])        
    elif set(dst_dims) == set(['band', 'day_obs', 'detector', 'physical_filter',
                             'tract', 'visit']):
        Dimension.update([
            ('visit', directories[-2]),
            ('physical_filter', directories[-3]),
            ('band', directories[-4]),
            ('day_obs', directories[-5]),
            ('tract', directories[-6]),
            ('detector',detectorIdComCam['_'.join(filenamedata.split('_')[7:9])])
        ])
    elif set(dst_dims) == set(['band', 'day_obs', 'detector',
                               'physical_filter', 'exposure']):
        Dimension.update([
            ('exposure', directories[-2]),
            ('day_obs', int(directories[-3])),
            ('band', filename.split('_')[3:4][0]),
            ('detector', detectorIdComCam['_'.join(filenamedata.split('_')[8:10])]),
            ('physical_filter', '_'.join(filename.split('_')[3:5]))
        ])
    elif set(dst_dims) == set(['visit', 'physical_filter', 'band',
                               'day_obs', 'detector']):
        Dimension.update([
            ('visit', int(directories[-2])),
            ('physical_filter', directories[-3]),
            ('band', directories[-4]),
            ('day_obs', int(directories[-5])),
            ('detector', detectorIdComCam['_'.join(filenamedata.split('_')[6:8])])
        ])
    elif set(dst_dims) == set(['visit', 'physical_filter', 'band', 'day_obs']):
        Dimension.update([
            ('visit', int(directories[-2])),
            ('physical_filter', directories[-3]),
            ('band', directories[-4]),
            ('day_obs', int(directories[-5]))
        ])
    elif set(dst_dims) == set(['band', 'patch', 'tract']):
        Dimension.update([
            ('band', directories[-2]),
            ('patch', int(directories[-3])),
            ('tract', int(directories[-4]))
        ])
    elif set(dst_dims) == set(['band', 'detector', 'physical_filter']):
        Dimension.update([
            ('physical_filter', directories[-2]),
            ('band', directories[-3]),
            ('detector', detectorIdComCam['_'.join(filenamedata.split('_')[5:7])])
        ])
    elif set(dst_dims) == set(['band', 'physical_filter', 'tract']):
        Dimension.update([
            ('physical_filter', directories[-2]),
            ('band', directories[-3]),
            ('tract', int(directories[-4]))
        ])
    elif set(dst_dims) == set(['patch', 'tract']):
        Dimension.update([
            ('patch', int(directories[-2])),
            ('tract', int(directories[-3]))
        ])
    elif set(dst_dims) == set(['band', 'tract']):
        Dimension.update([
            ('band', directories[-2]),
            ('tract', int(directories[-3]))
        ])
    elif set(dst_dims) == set(['tract']):
        Dimension.update([
            ('tract', int(directories[-2]))
        ])
    elif set(dst_dims) == set(['htm7']):
        Dimension.update([
            ('htm7', int(filename.split('.')[0]))
        ])
    elif set(dst_dims) == set(['band']):
        Dimension.update([
            ('band', directories[-2])
        ])
    elif set(dst_dims) == set(['detector']):
        Dimension.update([
            ('detector', detectorIdComCam['_'.join(filenamedata.split('_')[2:4])])
        ])
    elif dst_dims == [''] or dst_dims == []:
        pass
    else:
        print(f'[ERROR] dataType {dtype} not supported.')
        raise SystemExit
    if verbose > 3:
        print("<>&&<>",Dimension)
    return Dimension

#------------------------------------------------------------------------------

def create_csvs(dataDir, dataset_types):
    globPath = os.path.join(baseDir, dataDir, '**/*.*')
    print(f'Globbing {globPath} ...')
    filepaths = glob.glob(globPath, recursive=True)

    for dataType in dataset_types:
        print(f'Writing csv file for {dataType}...')
        csvPath = os.path.join(csvBaseDir, f'{dataType}_urls.csv')
        with open(csvPath, 'w') as f:
            for path in filepaths:
                if '/%s/' % dataType in path \
                   and not any(x in path for x in ['fallbackFlats']):
                    f.write(f'{path}\n')
        print(f'File list written to: {csvPath}')
    print('All dataType csv file lists written!')

#-----------------------------------------------------------------------------

def create_skymap(theSkymapFile):
    with open(skyMapPickleFile, 'r', errors='replace') as f:
        skymaplines = f.readlines()[1:]

    header = ['import lsst.skymap.ringsSkyMap\n',
              'import lsst.skymap.tractBuilder\n\n',
              "config.name='lsst_cells_v1'\n",
              "config.skyMap='rings'\n\n"]

    writeline = False    
    with open(theSkymapFile, 'w') as f:
        f.write(header[0])
        f.write(f'#{skymaplines[0]}')
        f.writelines(header[1:])
        for line in skymaplines[1:-1]:
            if writeline:
                f.write(line.replace("config.", "config.skyMap['rings']."))
            if line.startswith('import'):
                writeline = True 
    print(f'SkyMap written to: {theSkymapFile}')

#-----------------------------------------------------------------------------

def handle_exception(exception):
    print('An exception occurred:', str(exception))
    with open(f'autoIngest_error.log', 'a') as f:
        with redirect_stdout(f):
            traceback.print_exc()
            if "<WARNING>" in repr(exception): 
                f.write(f'{exception}\n')

#-----------------------------------------------------------------------------
        
def onRun(tasks, dataTypes=[], dataDir=defDataDir, testRun=False):
    """
    Run given tasks.
    """

    # create skymap config file
    if tasks['createSkymap']:
        create_skymap(skyMapConfigFile)

    # get all dataSetTypes
    allDataSetTypes = read_all_datatypes(dataTypesFile)

    if dataTypes:
        dataSetTypes = {k: allDataSetTypes[k]
                        for k in allDataSetTypes.keys() & set(dataTypes)}
    else:
        if tasks['Fromhere']:
            adstlist = list(allDataSetTypes.keys())
            adstsublist = adstlist[adstlist.index(tasks['Fromhere']):]
            dataSetTypes = {k: allDataSetTypes[k] for k in adstsublist} 
        else:
            dataSetTypes = allDataSetTypes

    # split dataset types into first and second ingest lists
    firstDataSetTypes = defaultdict(lambda: defaultdict(list))
    secondDataSetTypes = defaultdict(lambda: defaultdict(list))
    for k in dataSetTypes:
        if k in firstIngDsTypes[0]:
            firstDataSetTypes['coadd'][k] = dataSetTypes[k]
        elif k in firstIngDsTypes[1]:
            firstDataSetTypes['image'][k] = dataSetTypes[k]
        else:
            if not k in ['raw','skyMap', 'the_monster_20250219']:
                secondDataSetTypes[k] = dataSetTypes[k]
            
    # write csv files
    if tasks['makeList']:
        create_csvs(dataDir, dataSetTypes)

    # process each dataSetType
    if tasks['writeEcsv']:
        for dataType in dataSetTypes:
            write_ecsvs(dataType, dataSetTypes, dataDir, testRun)

    # setup the butler instrument and skymap
    if tasks['setupButler']:
        butler_setup(testRun)

    # ingest raw data and register visits
    if tasks['ingestraW']:
        butler_ingest_raw(testRun)

    # register coadd dataset types
    if tasks['Registerdata'] == 'c':
        butler_register_datasets(firstDataSetTypes['coadd'], testRun)

    # ingest coadd dataset types
    if tasks['Ingestdata'] == 'c':
        butler_ingest_datasets(firstDataSetTypes['coadd'], testRun)

    # define visits
    if tasks['defineVisits']:
        butler_define_visits(testRun)

    # register image dataset types
    if tasks['Registerdata'] == 'i':
        butler_register_datasets(firstDataSetTypes['image'], testRun)

    # ingest image dataset types
    if tasks['Ingestdata'] == 'i':
        butler_ingest_datasets(firstDataSetTypes['image'], testRun)

    # register all other dataset types
    if tasks['Registerdata'] == 'a':
        butler_register_datasets(secondDataSetTypes, testRun)

    # ingest all other dataset types
    if tasks['Ingestdata'] == 'a':
        butler_ingest_datasets(secondDataSetTypes, testRun)

    # register and ingest data one-by-one
    if tasks['regingdAta']:
        #butler_reg_ing_datasets(dataSetTypes, testRun)
        butler_reg_ing_datasets(secondDataSetTypes, testRun)
    # ingest calib
    if tasks['ingestCalib']:
        butler_ingest_calib(testRun)

    # ingest refcats
    if tasks['ingestMonster']:
        butler_ingest_refcats(testRun)

    # butler query database for already ingested data
    if tasks['Querybutler']:
        butler_query(tasks['Querybutler'], dataSetTypes, testRun)
        
#-----------------------------------------------------------------------------

def read_all_datatypes(theDataTypesFile):
    dataSetTypes = defaultdict(list)
    with open(theDataTypesFile, 'r') as f:
        datatypes = [line.strip().replace(' ','') for line in f.readlines()]

    for line in datatypes:
        theline = line.replace('DatasetType','').replace("'",'')
        dstelem = re.split(r',(?![^{}]*\})', theline[1:-1])
        dst, dims, format = dstelem[:3]
        dimensions = dims[1:-1].split(',')
        dataSetTypes[dst] = [dimensions, format]

    return dataSetTypes

#-----------------------------------------------------------------------------

def zrun_shell_cmd(cmdargs):
    # Start subprocess
    # bufsize = 1 means output is line buffered
    # universal_newlines = True is required for line buffering
    process = subprocess.Popen(cmdargs, shell=True,
                               bufsize=1,
                               stdout=PIPE,
                               stderr=PIPE,
                               universal_newlines=True)

    # Create callback function for process output
    buf = io.StringIO()
    def handle_output(stream, mask):
        # Because the process' output is line buffered, there's only ever one
        # line to read when this function is called
        line = stream.readline()
        buf.write(line)
        sys.stdout.flush()
        sys.stdout.write(line)

    # Register callback for an "available for read" event from subprocess'
    # stdout stream
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ, handle_output)

    # Loop until subprocess is terminated
    while process.poll() is None:
        # Wait for events and handle them with their registered callbacks
        events = selector.select()
        for key, mask in events:
            callback = key.data
            callback(key.fileobj, mask)

    # Ensure all remaining output is processed
    while True:
        line = process.stdout.readline()
        if not line:
            break
        buf.write(line)
        sys.stdout.write(line)


    errors = process.stderr.readlines()


    # Get process return code
    return_code = process.wait()
    selector.close()

    success = (return_code == 0)

    # Store buffered output
    output = buf.getvalue()
    buf.close()

    if errors:
        handle_exception(errors)
        raise SystemExit

    return output

#-----------------------------------------------------------------------------

def run_shell_cmd(cmdargs):
    """
    Run a system command.
    
    @param cmd:       The command to run.
    @type  cmd:       str

    @return: List of lines returned.
    @rtype:  list(str)

    """
    stdOut = []
    os.system(' '.join(cmdargs))
    #proc = subprocess.run(cmdargs, stdout=None, stderr=subprocess.PIPE)
    #proc = Popen(cmdargs, shell=True, stdout=STDOUT, stderr=PIPE, bufsize=1,
    #             universal_newlines=True)
    
    #for line in proc.stdout:
    #    print(line, end='')
    #    stdOut.append(line)
    #stdOut = 1
    #errors = proc.stderr.decode()#.readlines()
    #while True:
    #    data = proc.stdout.read(1024)   # Alternatively proc.stdout.read(1024)
    #    if len(data) == 0:
    #        break
    #    sys.stdout.flush()
    #    sys.stdout.write(data)
    #    stdOut.append(data)
    #

    #proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, bufsize=1)
    #stdOut = [x.rstrip() for x in proc.stdout]
    #if verbose > 10:
    #for x in stdOut:
    #for x in proc.stdout:
    #    print("[%s]> %s" % (str(mxTime.now()), x.rstrip()), flush=True)
    #    stdOut.append(x.rstrip())

    
    #if errors:
    #    handle_exception(errors)
    #    raise SystemExit
    #print(">>>",stdOut)
    return stdOut

#-----------------------------------------------------------------------------

def write_ecsvs(dataType, dataSetTypes, dataDir=defDataDir, testRun=False):

    # read dataType files into 'urls'
    csvPath = os.path.join(csvBaseDir, f'{dataType}_urls.csv')
    if verbose > 4:
        print("---",csvPath)
    
    urls = pd.read_csv(csvPath, header=None, names=['urls'])

    if urls.empty:
        handle_exception(f"<WARNING> no files in {csvPath}")
        return

    urls['filename'], urls['instrument'], urls['skymap'], urls['band'], \
    urls['patch'], urls['tract'], urls['visit'], urls['physical_filter'], \
    urls['day_obs'], urls['detector'], urls['htm7'], urls['exposure'] = zip(
        *urls['urls'].apply(lambda x: change_base_directory(x, dataType, dataSetTypes[dataType][0]).values()))

    useUrls = urls.copy(deep=True)

    if verbose > 4:
        print(dataType, "<>",dataSetTypes[dataType])
        
    # create groups
    # dataSetTypes[dst] = [dimensions, format]
    try:
        groups = ['filename']
        for dim in dataSetTypes[dataType][0]:
            if dim:
                groups.append(dim)

        if verbose > 4:     
            print(">>",groups)
    except KeyError:
        print("[ERROR] dataType doesn't exist in groups.")
        raise SystemExit

    outGroups = groups
    if verbose > 4:
        print("<<<",dataType,"::",outGroups)
        print(">>>",useUrls)
        
    # create one file containing all data
    new_urls = useUrls[outGroups]
    if verbose > 4:
        print(">>>",new_urls)

    file_table = Table.from_pandas(new_urls)
    if verbose > 2:
        print(file_table)

    # set output path for ecsv file
    ecsvPath = os.path.join(ecsvBaseDir, f'{dataType}.ecsv')

    if testRun:
        print(f'[TEST] write to single table: {ecsvPath}')
    else:
        print(f'Writing to single table: {ecsvPath}')
        file_table.write(ecsvPath, overwrite=True)

#-----------------------------------------------------------------------------

def usage(shopts, opts):
    cmdDesc = {
        's': 'create rings skymap config file from pickle file',
        'l': 'create dataset types csv file lists',
        'e': 'create ecsv files from csv files',
        'f': 'run "-e/-l" from this datasetType onwards',
        'b': 'setup the butler',
        'w': 'ingest the raw data',
        'r': 'register data: [c(oadds), i(mage), (a)ll others]',
        'i': 'ingest data: [c(oadds), i(mage), (a)ll others]',
        'a': 'register & ingest data one-by-one datasettype',
        'v': 'define the visits (run after coadds before images)',
        'c': 'ingest calib data',
        'm': 'ingest reference monster catalogues',
        'x': 'execute all commands in order',
        'h': 'print this help',
        'q': 'mode to butler query the database [i(nfo)/s(kymap/c(ollections)/t(ypes)/d(ata))]',
        'd': f'top level directory for the ingest files [default: {defDataDir}]',
        't': "test run, don't write/ingest data"
    }
    cmdOpts = []
    for s, o in zip(shopts, opts):
        cmdOpts.append(f'[-{s}/--{o}]')
    
    print(f"Usage: python autoIngest.py {' '.join(cmdOpts)} datasetTypes")
    print
    for s, o in zip(shopts, opts):
        print(f'    -{s}/--{o:18} {cmdDesc[s]}')
    print('  datasetTypes: the dataset types to be processed: optional, default use all')

#------------------------------------------------------------------------------

def main(argv):
    dataTypes = []
    testRun = False
    dataDir = defDataDir

    tasks = OrderedDict(
        [('createSkymap', False),
         ('makeList', False),
         ('writeEcsv', False),
         ('Fromhere', ''),
         ('setupButler', False),
         ('ingestraW', False),
         ('Registerdata', ''),
         ('Ingestdata', ''),
         ('regingdAta', False),
         ('defineVisits', False),
         ('ingestCalib', False),
         ('ingestMonster', False),
         ('Querybutler', ''),
         ('eXecuteall', False)])

    shortOpts = ''.join([c.lower() for s in tasks.keys()
                         for c in s if c.isupper()])
    longOpts = list(tasks.keys())
    shortOptsArg = ''
    for x in shortOpts:
        if x in ['q','f', 'r', 'i']:
            shortOptsArg += f'{x}:'
        else:
            shortOptsArg += x

    longOptsArg = []
    for lo in longOpts:
        if lo in ['Querybutler', 'Fromhere', 'Registerdata', 'Ingestdata']:
            longOptsArg.append(f'{lo}:')
        else:
            longOptsArg.append(lo)

    try:
        opts, args = getopt.getopt(
            argv[1:], 'hd:t' + shortOptsArg,
            ['help', 'datadir:', 'test'] + longOptsArg)

    except getopt.GetoptError:
        # print help information and exit:
        print(argv)
        usage('hdt' + shortOpts, ['help', 'datadir', 'test'] + longOpts) 
        raise SystemExit

    for o, a in opts:
        if o in ('-h','--help'):
            usage('hdt' + shortOpts, ['help', 'datadir', 'test'] + longOpts)
            raise SystemExit
        if o in ('-d', '--datadir'):
            dataDir = a
        if o in ('-t', '--test'):
            testRun = True
        if o in ('-s', '--createSkymap'):
            tasks['createSkymap'] = True
        if o in ('-l', '--makeList'):
            tasks['makeList'] = True
        if o in ('-e', '--writeEcsv'):
            tasks['writeEcsv'] = True
        if o in ('-f', '--fromHere'):
            tasks['Fromhere'] = a 
        if o in ('-b', '--setupButler'):
            tasks['setupButler'] = True
        if o in ('-w', '--ingestraW'):
            tasks['ingestraW'] = True
        if o in ('-r', '--Registerdata'):
            tasks['Registerdata'] = a
        if o in ('-i', '--Ingestdata'):
            tasks['Ingestdata'] = a
        if o in ('-a', '--regingdAta'):
            tasks["regingdAta"] = True
        if o in ('-v', '--defineVisits'):
            tasks['defineVisits'] = True
        if o in ('-c', '--ingestCalib'):
            tasks['ingestCalib'] = True
        if o in ('-m', '--ingestMonster'):
            tasks['ingestMonster'] = True
        if o in ('-q', '--Querybutler'):
            tasks['Querybutler'] = a.split(',')
        if o in ('-x', '--eXecuteall'):
            tasks['eXecuteall'] = True
            
    if len(args) == 0:
        dataTypes = []
    elif len(args) == 1:
        dataTypes = args[0].split(',')
    else:
        usage()
        raise SystemExit

    if tasks['eXecuteall']:
        for task in tasks:
            tasks[task] = True

    # create ecsv directory
    os.makedirs(ecsvBaseDir, exist_ok=True)

    try:
        onRun(tasks, dataTypes, dataDir, testRun)
    except Exception as e:
        handle_exception(e)

#------------------------------------------------------------------------------

if __name__ == '__main__':
    main(sys.argv)
